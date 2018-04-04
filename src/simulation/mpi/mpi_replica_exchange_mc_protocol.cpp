/*
 * mpi_replica_exchange_mc_protocol.cpp
 *
 *  Created on: 30 Sep 2013
 *      Author: jmacdona
 */
#include "mpi_replica_exchange_mc_protocol.h"


using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::pow;
using std::sqrt;
using std::asin;
using std::string;

using boost::lexical_cast;
using boost::split;
using boost::is_any_of;
using boost::trim;
using boost::optional;

namespace PRODART {
namespace POSE {
namespace SIM {
namespace MPI {


mpi_replica_exchange_mc_protocol_shared_ptr new_mpi_replica_exchange_mc_protocol( boost::mpi::communicator comm,
		unsigned long steps,
		MTRand::MTRand_shared_ptr mt_ptr){
	mpi_replica_exchange_mc_protocol_shared_ptr ptr(new mpi_replica_exchange_mc_protocol(comm, steps, 1.0, mt_ptr));
	return ptr;
}

bool mpi_replica_exchange_mc_protocol::load_replica_specs(std::string filename){


	typedef std::vector<std::string> string_vector;
	replica_specs.clear();
	replica_rank_map.clear();

	string lineStr;
	unsigned long length, lineNum = 0 ;
	string_vector SplitVec;
	//int pdbs_loaded = 0;

	ifstream in_file(filename.c_str(), ios::in);

	if (!in_file.is_open()){
		if (rep_comm.rank() == 0){
			cerr << "ERROR: mpi_replica_exchange_mc_protocol: could not load replica specs file: "
					<< filename << endl;
		}
		return false;
	}

	//cout << "reading file..." << endl;
	while ( !in_file.eof() ) {
		getline(in_file, lineStr);

		trim(lineStr);
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {

			split( SplitVec, lineStr, is_any_of("\t ") );

			string rep_potname = SplitVec[0];
			double rep_beta = lexical_cast<double>(SplitVec[1]);

			replica_spec_element ele;
			ele.potential_name = rep_potname;
			ele.beta = rep_beta;

			ele.pot_ptr = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container(rep_potname);

			if (!ele.pot_ptr){
				if (rep_comm.rank() == 0){
					cerr << "ERROR: mpi_replica_exchange_mc_protocol: problem loading potentials preset: "
							<< rep_potname
							<< endl;
				}
			}

			this->replica_specs.push_back(ele);

		}

	}

	if ((unsigned int)rep_comm.size() != (unsigned int)replica_specs.size()){
		if (rep_comm.rank() == 0){
			cerr << "ERROR: MPI communicator ("
					<< rep_comm.size()
					<< ") is not the same size as the replica specifications ("
					<< replica_specs.size()
					<< ")."
					<< " You probably need to rerun with corrected mpirun arguments."
					<< endl;
		}
		return false;
	}

	return true;
}

int mpi_replica_exchange_mc_protocol::get_num_replicas() const{
	return (int)replica_specs.size();
}





bool mpi_replica_exchange_mc_protocol::stage_initialise(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){

	state.replica_number = rep_comm.rank();

	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials =  replica_specs[state.replica_number].pot_ptr;
	protein->index();
	//protein->backup_coords();
	coords_backup = protein->get_all_coords();
	state.beta = this->beta;
	state.step_num = 0;
	state.was_big_move = false;
	state.stop = false;
	state.tries = 0;
	//meta_data->set_update_pair_lists_flag();
	state.steps_no_update = 0;
	meta_data->recalc_pair_lists_dists();
	state.energy_map = potentials->get_default_energies_map();
	state.energy = potentials->get_energy(meta_data, state.energy_map);
	state.prev_acc_energy = state.energy;
	state.best_energy = state.energy;
	state.start_energy_map = state.energy_map;
	state.best_energy_map = state.energy_map;
	state.start_energy = state.energy;
	state.best_pose = protein->clone();
	state.start_pose = protein->clone();
	state.cumulative_move = 0;
	potentials->init_set_up(meta_data, state.energy_map);

	if ((unsigned int)rep_comm.size() != (unsigned int)replica_specs.size()){
		if (rep_comm.rank() == 0){
			cerr << "ERROR: MPI communicator ("
					<< rep_comm.size()
					<< ") is not the same size as the replica specifications ("
					<< replica_specs.size()
					<< ")."
					<< " You probably need to rerun with corrected mpirun arguments."
					<< endl;
		}
		return false;
	}

	exch_swap_fwd_requests.clear();
	exch_swap_result_requests.clear();
	rec_swap_result_vals.clear();

	exch_swap_fwd_requests.resize((unsigned int)rep_comm.size());
	exch_swap_result_requests.resize((unsigned int)rep_comm.size());
	rec_swap_result_vals.resize((unsigned int)rep_comm.size());
	for (unsigned int i = 0; i < (unsigned int)rep_comm.size(); i++){
		exch_swap_fwd_requests[i] = request_shared_ptr();
		exch_swap_result_requests[i] = request_shared_ptr();
	}

	// only track rank_replica_map on rank 0
	if (rep_comm.rank() == 0){
		for (unsigned int i = 0; i < replica_specs.size(); i++){
			replica_rank_map[i] = i;
		}
	}
	//state.replica_number = rep_comm.rank();

	state.beta = replica_specs[state.replica_number].beta;

	exchange_msg_wait_flag = false;

	cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " initialised" << endl;

	return true;

}


bool mpi_replica_exchange_mc_protocol::stage_continue(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){


	if (state.stop == true){
		return false;
	}
	else if (state.step_num >= this->steps_to_run){
		return false;
	}
	else if (state.tries >= this->max_tries){
		return false;
	}

	return true;

}

bool mpi_replica_exchange_mc_protocol::stage_move(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials =  replica_specs[state.replica_number].pot_ptr;
	const bool result =  movers->make_move(meta_data);
	PRODART::POSE::MOVERS::mover_flags& flags = meta_data->get_mover_flags();
	state.was_big_move = flags.is_large_move;
	if (state.cumulative_move + flags.move_dist > this->max_move_no_update){
		state.was_big_move = true;
		//std::cout << "large move: " << state.step_num << "\t" << state.cumulative_move + flags.move_dist << "\n";
	}
	if (state.was_big_move || state.was_prev_big_move || state.steps_no_update >= this->max_steps_no_update) {
		//meta_data->set_update_pair_lists_flag();
		state.steps_no_update = 0;
		state.cumulative_move = 0;
	}
	meta_data->recalc_pair_lists_dists();
	state.energy = potentials->get_energy(meta_data, state.energy_map);




	return result;
}

bool mpi_replica_exchange_mc_protocol::stage_accept_test(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
    bool accept = false;
    const double deltaE = state.energy - state.prev_acc_energy;
    //std::cout << "deltaE:\t" << deltaE << std::endl;
    if (deltaE < 0) {
            accept = true;
    }
    else {
    	//std::cout << "beta:\t" << state.beta << std::endl;
    	const double prob_accept =  exp(-state.beta * deltaE);
    	//std::cout << "prob_accept:\t" << prob_accept << std::endl;
    	const double randnum = rand_gen->rand(1);
    	//std::cout << "randnum:\t" << randnum << std::endl;
    	if (prob_accept > randnum) accept = true;
    }



    //std::cout << "accept:\t" << accept << std::endl;
    return accept;
}

bool mpi_replica_exchange_mc_protocol::stage_accepted(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials =  replica_specs[state.replica_number].pot_ptr;
	state.prev_acc_energy = state.energy;
	state.prev_acc_energy_map = state.energy_map;
	PRODART::POSE::MOVERS::mover_flags& flags = meta_data->get_mover_flags();
	if (state.energy < state.best_energy){
		state.best_energy = state.energy;
		state.best_energy_map = state.energy_map;
		state.best_pose = protein->clone();
	}
	state.step_num++;
	state.steps_no_update++;
	//protein->backup_coords();
	coords_backup = protein->get_all_coords();
	state.tries = 0;
	state.was_big_move = false;
	state.was_prev_big_move = false;
	state.cumulative_move += flags.move_dist;
	potentials->accepted(meta_data, state.energy_map);

	/*
	std::cout << "accepted\t" << state.step_num << "\t";
	state.energy_map.print_weighted_components(std::cout);
	*/

	//protein->outputPdb(traj_out);

	return true;
}

bool mpi_replica_exchange_mc_protocol::stage_rejected(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials =  replica_specs[state.replica_number].pot_ptr;
	state.num_rej++;
	state.tries++;
	//state.steps_no_update++;
	//TODO investigate more efficient alternative to this for big moves:
	state.was_prev_big_move = state.was_big_move;
	//protein->restore_backed_up_coords();
	protein->restore_coords(coords_backup);
	if (state.tries > this->recalc_tries){
		//meta_data->set_update_pair_lists_flag();
		state.steps_no_update = 0;
		state.cumulative_move = 0;
		meta_data->recalc_pair_lists_dists();
		state.energy = potentials->get_energy(meta_data, state.energy_map);
	    const double deltaE = state.energy - state.prev_acc_energy;
		if (fabs(deltaE) > energy_tolerance){
			state.prev_acc_energy = state.energy;
			std::cerr << "\nrank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")"
					<< " ERROR: mpi_replica_exchange_mc_protocol: energy calculation problem\tdeltaE:\t" << deltaE << "\t"
					<< "at step:\t" << state.step_num
					<< std::endl;
			state.energy_map.print_unweighted_diff(cout, state.prev_acc_energy_map);
			//POSE::POTENTIALS::potentials_store diff_map =  state.energy_map.get_diff(state.prev_acc_energy_map);


			PRODART::UTILS::vector3d highest, lowest;
			PRODART::POSE_UTILS::get_ca_bounding_box(protein, lowest, highest);
			std::cout << "lowest:\t" << lowest << std::endl;
			std::cout << "highest:\t" << highest << std::endl;
		}
	}

	/*
	std::cout << "rejected\t" << state.step_num << "\t";
	state.energy_map.print_weighted_components(std::cout);
	*/

	return true;
}


bool mpi_replica_exchange_mc_protocol::stage_before_next_round(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){


	//MPI communication
	//cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " starting MPI communcations. Step: " << (state.num_rej + state.step_num) << endl;




	// receive stuff
	const int start_rank = (rep_comm.rank() == 0) ? 1 : 0;
	const int end_rank = (rep_comm.rank() == 0) ? rep_comm.size() - 1 : 0;

	/*
	cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " start_rank:" << start_rank << " end_rank:" << end_rank
			<< endl;
			*/
	for (int rank_ = start_rank; rank_ <= end_rank; rank_++){
		// post recieve request
		if (exch_swap_result_requests[rank_]){
			optional<boost::mpi::status> stat = exch_swap_result_requests[rank_]->test();
			if (stat){
				//recieved a message I think
				cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ") message received from rank:" << rank_
						<< " value: " << rec_swap_result_vals[rank_]
						<< endl;

				exch_receieve_type& rec_msg = rec_swap_result_vals[rank_];

				//! check if needs to be forwarded
				if (rep_comm.rank() == 0
						&& rec_msg.dest_replica_number != state.replica_number
						){
					//forward to destination replica
					{
						// check if dest_rank already set:
						const int dest_rank = rec_msg.dest_rank >= 0 ? rec_msg.dest_rank  : replica_rank_map[rec_msg.dest_replica_number];

						if (dest_rank != 0){
							rec_msg.dest_rank = dest_rank;
							cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " message received from rank:" << rank_
									<< " forwarding to rank:" <<  dest_rank
									<< endl;
							// should I keep track of request?
							exch_swap_fwd_requests[dest_rank] = request_shared_ptr(new boost::mpi::request);
							*(exch_swap_fwd_requests[dest_rank]) = rep_comm.isend(dest_rank, SWAP_SEND, rec_msg);
						}
						else if (rec_msg.confirm && (state.replica_number == rec_msg.source_replica_number)){
							cout << "hmmm should not be sending to myself but it is OK as a confirm message" << endl;
						}
						else {
							cout << "ERROR: hmmm should not be sending to myself" << endl;
						}

						if (rec_msg.confirm) {
							// TODO rank 0 housekeeping to swap replica rank mapping on confirmed swap
							cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")"
										<< " HOUSEKEEPING:me swap confirmed between replica:" << rec_msg.source_replica_number
										<< " at step " << rec_msg.source_step
										<< " and replica:" <<  rec_msg.dest_replica_number
										<< " at step " << rec_msg.dest_step
										<< endl;
							std::swap(replica_rank_map.find(rec_msg.source_replica_number)->second,
									replica_rank_map.find(rec_msg.dest_replica_number)->second);
							for (unsigned int i = 0; i < replica_rank_map.size(); i++ ){
								cout << i << "\t" << replica_rank_map[i] << endl;
							}
							if (dest_rank == 0){
								exchange_msg_confirm_flag = false;
							}
						}

					}
					//
				}
				else {



					if ((rec_msg.dest_replica_number == state.replica_number && !rec_msg.confirm)
							|| (rec_msg.source_replica_number == state.replica_number && rec_msg.confirm)){
						cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " replica:" << state.replica_number << " message reached destination - time to deal with it!"
								<< endl;

						//has reached destination replica (could be rank 0) so deal with it
						if (rec_msg.accepted && (!rec_msg.confirm) ){
							// accept message
							cout << "rank:" << rep_comm.rank() << " (replica:" << state.replica_number
									<< ") received acceptance message from replica:" <<  rec_msg.source_replica_number
									<< endl;
							// TODO clear waiting flags and do swap etc
							const int new_replica_num = rec_msg.source_replica_number;
							exchange_msg_wait_flag = false;
							{
								rec_msg.flip_and_confirm();
								exch_accept_request = request_shared_ptr(new boost::mpi::request);
								*(exch_accept_request) = rep_comm.isend(rank_, SWAP_SEND, rec_msg);

								if (rep_comm.rank() == 0) {
									// TODO rank 0 housekeeping to swap replica rank mapping on confirmed swap
									cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")"
													<< " HOUSEKEEPING:other swap confirmed between replica:" << rec_msg.source_replica_number
													<< " at step " << rec_msg.source_step
													<< " and replica:" <<  rec_msg.dest_replica_number
													<< " at step " << rec_msg.dest_step
													<< endl;
									std::swap(replica_rank_map.find(rec_msg.source_replica_number)->second,
											replica_rank_map.find(rec_msg.dest_replica_number)->second);
									for (unsigned int i = 0; i < replica_rank_map.size(); i++ ){
										cout << i << "\t" << replica_rank_map[i] << endl;
									}
								}
							}
							//
							state.replica_number = new_replica_num;//rec_msg.source_replica_number;
							state.beta = replica_specs[state.replica_number].beta;//;rec_msg.source_beta;
							state.step_num = replica_request_step;
							PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials =  replica_specs[state.replica_number].pot_ptr;
							protein->restore_coords(replica_request_coords_backup);
							meta_data->recalc_pair_lists_dists();
							state.energy_map = potentials->get_default_energies_map();
							state.prev_acc_energy = potentials->get_energy(meta_data, state.energy_map);
							state.prev_acc_energy_map = state.energy_map;
							cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << "  I am now replica:" << state.replica_number
									<< endl;
						}
						else if (rec_msg.rejected && (!rec_msg.confirm) ){
							// reject message
							cout << "rank:" << rep_comm.rank() << " (replica:" << state.replica_number
									<< ") received rejection message from replica:" <<  rec_msg.source_replica_number
									<< endl;
							// TODO clear waiting flags etc
							exchange_msg_wait_flag = false;
						}
						else if ((rec_msg.confirm) ){
							// accept or confirm message
							cout << "rank:" << rep_comm.rank() << " (replica:" << state.replica_number
									<< ") received confirmation message from replica:" <<  rec_msg.dest_replica_number
									<< endl;
							exchange_msg_confirm_flag = false;
						}
						else if ((!rec_msg.accepted) && (!rec_msg.confirm) ){

							const int other_replica = rec_msg.source_replica_number;

							if (exchange_msg_wait_flag || exchange_msg_confirm_flag){
								// waiting for my own swap request so reject
								cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " Rejected exchange request replica:" << state.replica_number
										<< " with replica:" << other_replica
										<< " because I am already waiting for my own request or confirmation"
										<< endl;
								rec_msg.flip_and_reject();
								exch_accept_request = request_shared_ptr(new boost::mpi::request);
								*(exch_accept_request) = rep_comm.isend(rank_, SWAP_SEND, rec_msg);
							}
							else {
								// metropolis criterion test
								double other_energy = state.prev_acc_energy;
								const bool pot_same = replica_specs[state.replica_number].potential_name.compare(replica_specs[other_replica].potential_name) == 0;
								const double beta_tol = 0.0000000000000001;
								const bool beta_same = fabs(replica_specs[state.replica_number].beta - rec_msg.source_beta) < beta_tol;
								if (pot_same == false){
									cout << "calculating energy" << endl;
									meta_data->recalc_pair_lists_dists();
									PRODART::POSE::POTENTIALS::potentials_energies_map temp_map = replica_specs[other_replica].pot_ptr->get_default_energies_map();
									other_energy = replica_specs[other_replica].pot_ptr->get_energy(meta_data, temp_map);
								}
								/* now apply metropolis criterion */
								if ((!pot_same) && beta_same){
									//Hamiltonian Replica Exchange
									// transition: (X_i, X'_j -> X'_i, X_j)
									// where X_i is structure X with potential U_i
									// and X'_j is structure X'_j with potential U_j
									// transition probability W(X_i, X'_j -> X'_i, X_j) = min(1,exp(-Delta))
									// Delta = beta * ( U_i(X') - U_j(X') + U_j(X) - U_i(X) )
									// X is this structure U_i is this potential
									// X' is other structure U_j is other potential
									//
									// U_i(X')  is this potential and other structure: rec_msg.dest_energy
									// U_j(X') is other potential and other structure: rec_msg.source_energy
									// U_j(X) is other potential and this structure: other_energy
									// U_i(X) is this potential and this structure:  state.prev_acc_energy
									const double delta = state.beta * (rec_msg.dest_energy - rec_msg.source_energy
											+ other_energy - state.prev_acc_energy);
									bool accept = false;
									double prob = 1;
									if (delta <= 0){
										accept = true;
									}
									else {
										prob = exp(-delta);
										const double randnum = rand_gen->rand(1);
										if (prob > randnum) accept = true;
									}

									if (accept){
										cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " Accepted HREM exchange replica:" << state.replica_number
												<< " with replica:" << other_replica
												<< " prob:" << prob
												<< endl;
										rec_msg.flip_and_accept();
										exchange_msg_confirm_flag = true;
										rec_msg.source_step = state.step_num;
										exch_accept_request = request_shared_ptr(new boost::mpi::request);
										*(exch_accept_request) = rep_comm.isend(rank_, SWAP_SEND, rec_msg);
										// TODO do swap

										state.replica_number = other_replica;
										state.beta = replica_specs[state.replica_number].beta;//rec_msg.dest_beta;
										PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials =  replica_specs[state.replica_number].pot_ptr;
										meta_data->recalc_pair_lists_dists();
										state.energy_map = potentials->get_default_energies_map();
										state.prev_acc_energy = potentials->get_energy(meta_data, state.energy_map);
										state.prev_acc_energy_map = state.energy_map;

										cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << "  I am now replica:" << state.replica_number
												<< endl;
									}
									else {
										cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " Rejected HREM exchange replica:" << state.replica_number
												<< " with replica:" << other_replica
												<< " prob:" << prob
												<< endl;
										rec_msg.flip_and_reject();
										exch_accept_request = request_shared_ptr(new boost::mpi::request);
										*(exch_accept_request) = rep_comm.isend(rank_, SWAP_SEND, rec_msg);

									}
								}
								else if (pot_same){
									//Temperature Replica Exchange
									// transition (X beta_X, X' beta_X' -> X beta_X', X' beta_X)
									// W(X beta_X, X' beta_X' -> X beta_X', X' beta_X) = min(1, exp(-Delta))
									// Delta = (beta_X - beta_X')(U(X') - U(X))
									const double delta = (state.beta - rec_msg.source_beta)
														* (rec_msg.source_energy - state.prev_acc_energy);

									bool accept = false;
									double prob = 1;
									if (delta <= 0){
										accept = true;
									}
									else {
										prob = exp(-delta);
										const double randnum = rand_gen->rand(1);
										if (prob > randnum) accept = true;
									}

									if (accept){
										cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " Accepted TREM exchange replica:" << state.replica_number
												<< " energy: " << state.prev_acc_energy
												<< " with replica:" << other_replica
												<< " energy " << rec_msg.source_energy
												<< " prob: " << prob
												<< " delta: " << delta
												<< " state.beta: " << state.beta
												<< " rec_msg.source_beta " << rec_msg.source_beta
												<< endl;
										rec_msg.flip_and_accept();
										exchange_msg_confirm_flag = true;
										rec_msg.source_step = state.step_num;
										exch_accept_request = request_shared_ptr(new boost::mpi::request);
										*(exch_accept_request) = rep_comm.isend(rank_, SWAP_SEND, rec_msg);
										// TODO do swap

										state.replica_number = other_replica;
										state.beta = rec_msg.dest_beta;
										cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << "  I am now replica:" << state.replica_number
												<< endl;
									}
									else {
										cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " Rejected TREM exchange replica:" << state.replica_number
												<< " with replica:" << other_replica
												<< " prob:" << prob
												<< endl;
										rec_msg.flip_and_reject();
										exch_accept_request = request_shared_ptr(new boost::mpi::request);
										*(exch_accept_request) = rep_comm.isend(rank_, SWAP_SEND, rec_msg);
									}

								}
								else {
									cout << "ERROR: unsupported Temperature Hamiltonian Replica Exchange"
											<< " " << replica_specs[state.replica_number].potential_name
											<< " " << replica_specs[other_replica].potential_name
											<< " " << replica_specs[state.replica_number].beta
											<< " " << replica_specs[other_replica].beta
											<< endl;
								}


								///////
							}
						}
					}
					else {
						// message received by wrong replica or another exchange happened in meantime
						/****************************************************************************
						 *
						 ****************************************************************************
						 */
						cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " WARNING: may have reached wrong destination (I have just changed replica_number)"
								<< endl;
						if (rec_msg.accepted == false && rec_msg.confirm == false){
							rec_msg.flip_and_reject();
							exch_accept_request = request_shared_ptr(new boost::mpi::request);
							*(exch_accept_request) = rep_comm.isend(rank_, SWAP_SEND, rec_msg);
							cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " sent rejection" << endl;
						}
						else {
							cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " ERROR!!!!!!!!!!!!!!!!! Reached wrong destination and I don't know how to deal with it"
								<< endl;
						}
					}
				}

				exch_swap_result_requests[rank_] = request_shared_ptr();
				rec_msg = exch_receieve_type();
			}
			else {
				//cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " no message got " << rec_msg << endl;
				//rec_msg = exch_receieve_type();
			}
		}
		//cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " passed 1" << endl;
		if (!exch_swap_result_requests[rank_]){
			// post receive request
			cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " posting receive request from rank:" << rank_
					<< " of type " << typeid(rec_swap_result_vals[rank_]).name()
					<< endl;
			exch_swap_result_requests[rank_] = request_shared_ptr(new boost::mpi::request);
			*(exch_swap_result_requests[rank_]) = rep_comm.irecv(rank_, SWAP_SEND, rec_swap_result_vals[rank_]);
			//exch_swap_result_requests[rank_]->wait();
			/*
			cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " waited for message received from rank:" << rank_
					<< " value: " << rec_swap_result_vals[rank_]
					<< endl;
			 */

		}
	}


	//cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " passed 2" << endl;

	// send stuff
	//bool bypass = false;//true;
	if (((state.num_rej + state.step_num) % this->exchange_frequency == 0)){

		if (rep_comm.rank() == 0){
			for (unsigned int i = 0; i < replica_rank_map.size(); i++ ){
				cout << "CENTRAL_CHECK " << i << "\t" << replica_rank_map[i] << endl;
			}
		}
		else {
			cout << "LOCAL_CHECK " << state.replica_number << "\t" << rep_comm.rank() << endl;
		}

		if (exch_swap_request){
			optional<boost::mpi::status> stat = exch_swap_request->test();
			if (stat){
				exch_swap_request = request_shared_ptr();
				cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " successful send" << endl;
			}
			else {
				cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " UNsuccessful send" << endl;
				// TODO NOTE not sure if I should do this.
				exch_swap_request->cancel();
				exch_swap_request = request_shared_ptr();
				//exchange_msg_confirm_flag = false;
			}
		}

		if (!exch_swap_request){
			if ((state.replica_number+1) < rep_comm.size()){
				const int next_replica_up = (state.replica_number+1) < rep_comm.size() ? state.replica_number+1 : 0;
				cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " replica:" << state.replica_number << " thinking about exchanging at step " <<  state.step_num
						<< " with replica:" << next_replica_up
						<< endl;
				if (exchange_msg_wait_flag){
					cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " step " <<  state.step_num
							<< " but wait flag is set"
							<< endl;
				}
				if (exchange_msg_confirm_flag){
					cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " step " <<  state.step_num
							<< " but confirm flag is set"
							<< endl;
				}
				if ((!exchange_msg_wait_flag) && (!exchange_msg_confirm_flag)){
					double dest_beta = replica_specs[next_replica_up].beta;

					double dest_energy = state.prev_acc_energy;
					if (replica_specs[state.replica_number].potential_name.compare(replica_specs[next_replica_up].potential_name) != 0){
						cout << "calculating energy" << endl;
						meta_data->recalc_pair_lists_dists();
						PRODART::POSE::POTENTIALS::potentials_energies_map temp_map = replica_specs[next_replica_up].pot_ptr->get_default_energies_map();
						dest_energy = replica_specs[next_replica_up].pot_ptr->get_energy(meta_data, temp_map);
					}
					if (rep_comm.rank() != 0){
						const int dest_rank = 0 ;
						exch_receieve_type msg( state.replica_number, //source_replica_number_,
								next_replica_up, //dest_replica_number_,
								state.beta, //double source_beta_,
								state.prev_acc_energy, //double source_energy_,
								dest_beta,//double dest_beta_,
								dest_energy,
								state.step_num,
								rep_comm.rank());//double dest_energy_);
						cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " sending message to rank:" << dest_rank << " of type: " << typeid(msg).name() << endl;
						exch_swap_request = request_shared_ptr(new boost::mpi::request);
						*(exch_swap_request) = rep_comm.isend(dest_rank, SWAP_SEND, msg);
						cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " sent" << endl;
					}
					else {
						exch_receieve_type msg( state.replica_number, //source_replica_number_,
								next_replica_up, //dest_replica_number_,
								state.beta, //double source_beta_,
								state.prev_acc_energy, //double source_energy_,
								dest_beta,//double dest_beta_,
								dest_energy,
								state.step_num,
								rep_comm.rank());//double dest_energy_);
						const int dest_rank = replica_rank_map[next_replica_up];

						cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " sending message to rank:" << dest_rank << " of type: " << typeid(msg).name() << endl;
						exch_swap_request = request_shared_ptr(new boost::mpi::request);
						*(exch_swap_request) = rep_comm.isend(dest_rank, SWAP_SEND, msg);
						cout << "rank:" << rep_comm.rank() << "(replica:" << state.replica_number << ")" << " sent" << endl;
					}
					replica_request_coords_backup = protein->get_all_coords();
					replica_request_energy = state.prev_acc_energy;
					replica_request_step = state.step_num;
					exchange_msg_wait_flag = true;

				}
				//
			}
		}
		else {
			//pending send
		}
	}


	return true;
}




bool mpi_replica_exchange_mc_protocol::stage_finish(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
	state.final_pose = protein->clone();
	//rep_comm.barrier();
	return true;
}







}
}
}
}

