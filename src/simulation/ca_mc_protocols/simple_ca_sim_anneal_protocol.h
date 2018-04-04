/*
 * simple_ca_sim_anneal_protocol.h
 *
 *  Created on: 15 Oct 2010
 *      Author: jmacdona
 */

#ifndef SIMPLE_CA_SIM_ANNEAL_PROTOCOL_H_
#define SIMPLE_CA_SIM_ANNEAL_PROTOCOL_H_
#include "simulation/mc_protocol_interface.h"
#include <limits>

namespace PRODART {
namespace POSE {
namespace SIM {
namespace CA {


class simple_ca_sim_anneal_protocol;


typedef boost::shared_ptr<simple_ca_sim_anneal_protocol> simple_ca_sim_anneal_protocol_shared_ptr;

mc_protocol_interface_shared_ptr new_simple_ca_sim_anneal_protocol(unsigned long start_beta_steps_,
		unsigned long anneal_steps_,
		unsigned long final_beta_steps_,
		double start_beta_,
		double final_beta_,
		MTRand::MTRand_shared_ptr mt_ptr);

class simple_ca_sim_anneal_protocol : public mc_protocol_interface {


	friend mc_protocol_interface_shared_ptr new_simple_ca_sim_anneal_protocol(unsigned long start_beta_steps_,
			unsigned long anneal_steps_,
			unsigned long final_beta_steps_,
			double start_beta_,
			double final_beta_,
			MTRand::MTRand_shared_ptr mt_ptr);

private:


protected:


	unsigned long start_beta_steps, anneal_steps, final_beta_steps;
	double start_beta, final_beta;

	simple_ca_sim_anneal_protocol(unsigned long start_beta_steps_,
			unsigned long anneal_steps_,
			unsigned long final_beta_steps_,
			double start_beta_,
			double final_beta_,
			MTRand::MTRand_shared_ptr mt_ptr) : mc_protocol_interface(start_beta_steps_ + anneal_steps_ + final_beta_steps_, start_beta_, mt_ptr),
			start_beta_steps(start_beta_steps_),
			anneal_steps(anneal_steps_),
			final_beta_steps(final_beta_steps_),
			start_beta(start_beta_),
			final_beta(final_beta_) {

	}

public:

	inline bool stage_accepted(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){
		state.prev_acc_energy = state.energy;
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


		//std::cout << "accepted\t" << state.step_num << "\t" << state.beta << "\t";
		//state.energy_map.print_weighted_components(std::cout);


		//protein->outputPdb(traj_out);

		if (state.step_num < this->start_beta_steps ){
			state.beta = this->start_beta;
		}
		else if (state.step_num >= this->start_beta_steps
				&& state.step_num < start_beta_steps + anneal_steps){
			//increment beta
			const double incr = (final_beta - start_beta) / (double)anneal_steps;
			state.beta = this->start_beta + (incr * (double)(state.step_num - this->start_beta_steps));
		}
		else if (state.step_num >= start_beta_steps + anneal_steps ){
			state.beta = this->final_beta;
		}

		return true;
	}


};





}
}
}
}

#endif /* SIMPLE_CA_SIM_ANNEAL_PROTOCOL_H_ */
