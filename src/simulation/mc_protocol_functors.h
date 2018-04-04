/*
 * mc_protocol_functors.h
 *
 *  Created on: 16 Sep 2013
 *      Author: jmacdona
 */

#ifndef MC_PROTOCOL_FUNCTORS_H_
#define MC_PROTOCOL_FUNCTORS_H_

#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include "pose_meta/pose_meta_interface.h"
#include "pose_meta/ca_pose_meta.h"
#include "potentials/potentials_name.h"
#include "potentials/potentials_store.h"
#include "potentials/potential_interface.h"
#include "potentials/potentials_container.h"
#include "prodart_env/prodart_env.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"


namespace PRODART {
namespace POSE {
namespace SIM {





class default_mc_state{

private:

protected:

public:

	unsigned long step_num;
	unsigned long steps_no_update;
	double cumulative_move;
	double start_energy;
	double energy;
	double best_energy;
	PRODART::POSE::POTENTIALS::potentials_energies_map start_energy_map;
	PRODART::POSE::POTENTIALS::potentials_energies_map energy_map;
	PRODART::POSE::POTENTIALS::potentials_energies_map prev_acc_energy_map;
	PRODART::POSE::POTENTIALS::potentials_energies_map best_energy_map;
	PRODART::POSE::pose_shared_ptr start_pose;
	PRODART::POSE::pose_shared_ptr final_pose;
	PRODART::POSE::pose_shared_ptr best_pose;
	double prev_acc_energy;
	unsigned long num_rej;
	unsigned long tries;
	double beta; // thermo beta
	bool stop;
	bool was_big_move;
	bool was_prev_big_move;
	int replica_number;



	default_mc_state() : step_num(0),
			steps_no_update(0),
			cumulative_move(0),
			start_energy(std::numeric_limits<double>::max()),
			energy(std::numeric_limits<double>::max()),
			best_energy(std::numeric_limits<double>::max()),
			start_pose(),
			final_pose(),
			best_pose(),
			prev_acc_energy(std::numeric_limits<double>::max()),
			num_rej(0),
			tries(0),
			beta(1.0),
			stop(false),
			was_big_move(false),
			was_prev_big_move(false),
			replica_number(-1){

	}

	std::ostream& print_summary(std::ostream& output);


};

inline std::ostream& default_mc_state::print_summary(std::ostream& output){
	output << "\nMC Summary:\n";
	if (replica_number >= 0){
		output << "\nreplica_number:\t" << this->replica_number << "\n";
	}
	output << "\nsteps:\t" << this->step_num << "\n";
	output << "beta:\t" << this->beta << "\n";
	output << "acceptance_rate:\t" << (static_cast<double>(step_num) / static_cast<double>(step_num + num_rej))   << "\n";
	output << "\n";
	output << "potential_name:\t"; this->start_energy_map.print_headers(output);
	output << "weights:\t"; this->start_energy_map.print_weights(output);
	output << "start:\t"; this->start_energy_map.print_weighted_components(output);
	output << "end:\t"; this->energy_map.print_weighted_components(output);
	output << "best:\t"; this->best_energy_map.print_weighted_components(output);
	output << "\n";
	if (this->best_pose && this->start_pose && this->final_pose){
		output << "CA_RMSDs:\tstart-end\tstart-best\tend-best\n";
		output << "CA_RMSDs:\t" << PRODART::POSE_UTILS::get_ca_rmsd(this->start_pose, this->final_pose) << "\t"
								<< PRODART::POSE_UTILS::get_ca_rmsd(this->start_pose, this->best_pose) << "\t"
								<< PRODART::POSE_UTILS::get_ca_rmsd(this->final_pose, this->best_pose)
								<< "\n";
		output << "\n";
		output << "CA_RadGyr:\tchain\tstart\tend\tbest\n" ;
		const int num_chains = this->start_pose->get_chain_count();
		for (int i = 0; i<num_chains;i++){
			output << "CA_RadGyr:\t" << this->start_pose->get_chain(i)->getChainID() << "\t"
					<< PRODART::POSE_UTILS::get_ca_radgyr(this->start_pose, i) << "\t"
					<< PRODART::POSE_UTILS::get_ca_radgyr(this->final_pose, i) << "\t"
					<< PRODART::POSE_UTILS::get_ca_radgyr(this->best_pose, i)
					<< "\n";
		}
	}
	output << std::endl;
	return output;
}






}
}
}
#endif /* MC_PROTOCOL_FUNCTORS_H_ */
