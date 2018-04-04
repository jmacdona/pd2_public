/*
 * structure_prediction.cpp
 *
 *  Created on: 9 Oct 2013
 *      Author: jmacdona
 */
#include "structure_prediction.h"
#include "protocols/protocols.h"

using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::SIM;

namespace PRODART {
namespace PROTOCOLS{
namespace MPI {
namespace STRUCTURE_PREDICTION {

//taxicab distance
bool has_neighbouring_contact(PRODART::POSE::META::ca_pose_meta_shared_ptr meta_,
		const int res1, const int res2,
		const int taxicab_dist){

	pose_shared_ptr pose_ = meta_->get_pose();
	const int rescount = pose_->get_residue_count();


	for (int i = res1-taxicab_dist; i <= res1+taxicab_dist;i++){
		for (int j = res2-taxicab_dist; j <= res2+taxicab_dist; j++){
			if (i >= 0 && i < rescount && j >=0 && j < rescount && i != j){
				const int dist = abs(i -res1) + abs(j - res2);
				if (dist <= taxicab_dist){
					upper_lower_bonded_pair_element ele = meta_->get_ca_GO_restraint(min(i,j),max(i,j));
					if (ele.is_valid == true && ele.half_bond_const > 0){
						return true;
					}

				}
			}
		}
	}


	return false;
}

bool mpi_replica_exchange_folding(PRODART::POSE::residue_type_vector seq,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const std::string psicov_file,
		const unsigned long REM_steps,
		boost::mpi::communicator comm,
		const std::string rep_specs,
		pose_shared_ptr start_struct,
		const unsigned long min_steps,
		const std::string bb_min_potentials_weight_set ){

	const int exchange_freq = 1000;
	const char chainID = 'A';
	const int start_res = 1;
	const int seq_sep_cutoff = 5;
	const double pred_con_energy = 1.0;
	const double rescon_half_bond_const = 1.0 / 64.0;
	const double non_native_energy = -0.1;
	const double psicov_prob_cutoff = 0.0;
	double secs_rst_wt = 5.0;

	pose_shared_ptr pose_;

	if (start_struct){
		cout << "INFO: starting from input structure" << endl;
		pose_ = start_struct;
		pose_->renumber_residues(1);
		pose_->get_chain(0)->setChainID('A');
	}
	else {
		cout << "INFO: starting from extended structure" << endl;
		pose_ = MISC::make_extended_ca_chain(seq, chainID, start_res, 0.2);
	}

	bool result = CA2MAIN::auto_add_secs_restraints(pose_, secs, secs_rst_wt);

	ifstream psicov_if(psicov_file.c_str(), ios::in);

	if (psicov_if.is_open()){
		std::map< int,std::map< int, boost::tuple<double, double, double> > > psicov_contacts = MISC::parse_psicov_file(psicov_if);
		typedef std::map< int,std::map< int, boost::tuple<double, double, double> > >::iterator outer_it;
		typedef std::map< int, boost::tuple<double, double, double> >::iterator inner_it;
		// apply psicov restraints
		//

		for (outer_it it_i = psicov_contacts.begin(); it_i != psicov_contacts.end(); it_i++){
			const int res_i = it_i->first;
			for (inner_it it_j = it_i->second.begin(); it_j != it_i->second.end(); it_j++){
				const int res_j = it_j->first;
				const int seq_sep = abs(res_j - res_i);
				if (it_j->second.get<2>() >= psicov_prob_cutoff && seq_sep >= seq_sep_cutoff){
					restraints_store::Instance()->add_ca_GO_rst(lexical_cast<string>(res_i), chainID,
							lexical_cast<string>(res_j), chainID,
							0,12, pred_con_energy);
					restraints_store::Instance()->add_contact_rst(lexical_cast<string>(res_i), chainID,
							lexical_cast<string>(res_j), chainID,
							0, 12, rescon_half_bond_const);
				}
				/* deal with this in later loop
				else {
					restraints_store::Instance()->add_ca_GO_rst(lexical_cast<string>(res_i), chainID,
							lexical_cast<string>(res_j), chainID,
							0,12, non_native_energy);
				}
				*/

			}
		}




		potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot = pot_factory->make_preset_potentials_container("ca_default");

		PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(pose_);
		//check restraints
		for (int i = 0; i < pose_->get_residue_count(); i++){
			for (int j = i+1; j < pose_->get_residue_count(); j++){
				upper_lower_bonded_pair_element ele = meta_data->get_ca_GO_restraint(i,j);
				if (j-i <= seq_sep_cutoff){
					ele.half_bond_const = 0;
					ele.is_valid = true;
					meta_data->add_ca_GO_restraint(i, j, ele);
				}

				else if ( ele.is_valid == false && has_neighbouring_contact(meta_data, i, j, 4) ){
					ele.half_bond_const = 0;
					ele.is_valid = true;
					meta_data->add_ca_GO_restraint(i, j, ele);
				}

				else if (ele.is_valid == false){
					ele.half_bond_const = non_native_energy;
					ele.equilibrium_dist_lower = 0.0;
					ele.equilibrium_dist_upper = 12.0;
					meta_data->add_ca_GO_restraint(i, j, ele);
					/*
					cout << "adding non-native contact energy " << i+1 << " "
							<< j+1 << " "
							<< non_native_energy
							<< endl;
					 */
				}
			}
		}


		PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = MOVERS::move_set_factory::Instance()->make_preset_move_set(meta_data,
				"ca_standard",
				PRODART::ENV::get_random_num_gen());

		PRODART::POSE::SIM::MPI::mpi_replica_exchange_mc_protocol_verbose_shared_ptr protocol =
				PRODART::POSE::SIM::MPI::new_mpi_replica_exchange_mc_protocol_verbose( comm,
				REM_steps,
				ENV::get_random_num_gen());

		protocol->set_exchange_frequency(exchange_freq);

		PRODART::POSE::SIM::monte_carlo_shared_ptr mc_sim = PRODART::POSE::SIM::new_monte_carlo(pot,
				overall_move_set,
				protocol);

		if (protocol->load_replica_specs(rep_specs)){

			cout << "running ca anneal" << endl;

			mc_sim->make_move(meta_data);

			const int replica_number = mc_sim->get_state().replica_number;

			string output_pdb_filename = PRODART::ENV::get_option_value<string>("output_root");
			output_pdb_filename.append(".replica");
			output_pdb_filename.append(boost::lexical_cast<string>(replica_number));
			output_pdb_filename.append(".pdb");

			ofstream output_pdb_file(output_pdb_filename.c_str(), ios::out);
			if (!output_pdb_file.is_open()){
				cerr << "ERROR: can't open output file: " << output_pdb_filename
						<< endl;
				return -1;
			}

			pose_->outputPdb(output_pdb_file, true);

			mc_sim->get_state().print_summary(cout);
		}
		else {
			cerr << "ERROR: mpi_replica_exchange_folding: some problem with replica specs file or mpirun arguments: " << ENV::get_option_value<string>("mpi_mc:replica_specs")
							<< endl;
		}
	}
	else {
		// can't open psicov input
		return false;
	}



	return true;
}








}
}
}
}
