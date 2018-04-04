/*
 * mpi_replica_exchange_mc_protocol_verbose.h
 *
 *  Created on: 3 Oct 2013
 *      Author: jmacdona
 */

#ifndef MPI_REPLICA_EXCHANGE_MC_PROTOCOL_VERBOSE_H_
#define MPI_REPLICA_EXCHANGE_MC_PROTOCOL_VERBOSE_H_



#include "mpi_replica_exchange_mc_protocol.h"


namespace PRODART {
namespace POSE {
namespace SIM {
namespace MPI {

class mpi_replica_exchange_mc_protocol_verbose;
typedef boost::shared_ptr<mpi_replica_exchange_mc_protocol_verbose> mpi_replica_exchange_mc_protocol_verbose_shared_ptr;
mpi_replica_exchange_mc_protocol_verbose_shared_ptr new_mpi_replica_exchange_mc_protocol_verbose( boost::mpi::communicator comm,
		unsigned long steps,
		MTRand::MTRand_shared_ptr mt_ptr);

class mpi_replica_exchange_mc_protocol_verbose : public mpi_replica_exchange_mc_protocol {

	friend mpi_replica_exchange_mc_protocol_verbose_shared_ptr new_mpi_replica_exchange_mc_protocol_verbose( boost::mpi::communicator comm,
			unsigned long steps,
			MTRand::MTRand_shared_ptr mt_ptr);

protected:
	mpi_replica_exchange_mc_protocol_verbose(boost::mpi::communicator comm,
			unsigned long steps,
			double _beta,
			MTRand::MTRand_shared_ptr mt_ptr) : mpi_replica_exchange_mc_protocol(comm,
					steps,
					_beta,
					mt_ptr),
					traj_out_freq(1000) {}


	unsigned long traj_out_freq;
	std::vector< boost::shared_ptr<std::ostream> > traj_out_streams;

public:

	void stage_initialise_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){

		traj_out_streams.clear();
		const int num_procs = rep_comm.size();
		for (int i = 0; i < num_procs; i++){
			std::string filename = PRODART::ENV::get_option_value<std::string>("output_root");
			filename.append("_traj");
			filename.append(".replica");
			filename.append(boost::lexical_cast<std::string>(i));
			filename.append(".pdb");
			{
				std::ofstream clearout(filename.c_str(), std::ios::out);
				if (clearout.is_open()){
					clearout.close();
				}
			}
			traj_out_streams.push_back(boost::shared_ptr<std::ostream>(new std::ofstream(filename.c_str(), std::ios::app)));
		}

	}

	void stage_accepted_print(default_mc_state& state,
				const PRODART::POSE::pose_shared_ptr protein,
				PRODART::POSE::META::pose_meta_shared_ptr meta_data,
				PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
				PRODART::POSE::MOVERS::mover_shared_ptr movers){

		if (!exchange_msg_wait_flag){
			if (state.step_num % traj_out_freq == 0){
				protein->outputPdb(*(traj_out_streams[state.replica_number]), true);

				std::cout << "rank:" << rep_comm.rank() << " replica:" << state.replica_number
						<< " step:" << state.step_num
						<< " components:\t";
				state.energy_map.print_headers(std::cout);

				std::cout << "rank:" << rep_comm.rank() << " replica:" << state.replica_number
						<< " step:" << state.step_num
						<< " ENERGY:\t";
				state.energy_map.print_weighted_components(std::cout);

			}
		}
		else {
			if (state.step_num % traj_out_freq == 0){

				std::cout << "rank:" << rep_comm.rank() << " replica:" << state.replica_number
						<< " skipping traj output as exchange_msg_wait_flag is set"
						<< std::endl;
			}
		}

	}

};



}
}
}
}

#endif /* MPI_REPLICA_EXCHANGE_MC_PROTOCOL_VERBOSE_H_ */
