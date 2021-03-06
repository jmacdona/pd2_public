/*
 * ca_motif_anneal_protocol.h
 *
 *  Created on: 5 Jan 2011
 *      Author: jmacdona
 */

#ifndef CA_MOTIF_ANNEAL_PROTOCOL_H_
#define CA_MOTIF_ANNEAL_PROTOCOL_H_

#include "simulation/mc_protocol_interface.h"
#include <limits>
#include <iostream>
#include <fstream>
#include "potentials/potentials_factory.h"

namespace PRODART {
namespace POSE {
namespace SIM {
namespace CA {


class ca_motif_anneal_protocol;


typedef boost::shared_ptr<ca_motif_anneal_protocol> ca_motif_anneal_protocol_shared_ptr;

mc_protocol_interface_shared_ptr new_ca_motif_anneal_protocol(unsigned long start_beta_steps_,
		unsigned long anneal_steps_,
		unsigned long final_beta_steps_,
		double start_beta_,
		double final_beta_,
		MTRand::MTRand_shared_ptr mt_ptr);

class ca_motif_anneal_protocol : public mc_protocol_interface {


	friend mc_protocol_interface_shared_ptr new_ca_motif_anneal_protocol(unsigned long start_beta_steps_,
			unsigned long anneal_steps_,
			unsigned long final_beta_steps_,
			double start_beta_,
			double final_beta_,
			MTRand::MTRand_shared_ptr mt_ptr);

private:

	//std::ofstream out_traj;

protected:


	unsigned long start_beta_steps, anneal_steps, final_beta_steps;
	double start_beta, final_beta;
	POTENTIALS::as_geom_pot_shared_ptr as_pot;
	double motif_move_prob;
	bool motif_loaded;
	bool debug;
	POTENTIALS::int_int_map_vector store_motif_mapping;

	ca_motif_anneal_protocol(unsigned long start_beta_steps_,
			unsigned long anneal_steps_,
			unsigned long final_beta_steps_,
			double start_beta_,
			double final_beta_,
			MTRand::MTRand_shared_ptr mt_ptr) : mc_protocol_interface(start_beta_steps_ + anneal_steps_ + final_beta_steps_, start_beta_, mt_ptr),
			start_beta_steps(start_beta_steps_),
			anneal_steps(anneal_steps_),
			final_beta_steps(final_beta_steps_),
			start_beta(start_beta_),
			final_beta(final_beta_){
		as_pot = boost::static_pointer_cast<POTENTIALS::as_geom_pot, POTENTIALS::potential_interface>(PRODART::POSE::POTENTIALS::potentials_factory::Instance()->get_potential(POTENTIALS::potentials_name("as_geom_pot")));
		motif_move_prob = 0.1;
		motif_loaded = false;
		debug = false;
	}

public:


	inline bool stage_initialise(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){

		motif_loaded = as_pot->is_loaded();
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
		//as_pot->random_assign_motif(meta_data->get_pose());
		if (motif_loaded){
			//as_pot->move_store_mapping();
			store_motif_mapping = as_pot->get_motif_assignment();
		}
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

		/*
		state.energy_map.print_headers(std::cout);
		state.energy_map.print_weights(std::cout);
		std::cout << "start\t" << state.step_num << "\t";
		state.energy_map.print_weighted_components(std::cout);
		*/

		return true;

	}

	inline bool stage_move(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){
		const bool result =  movers->make_move(meta_data);


		if (motif_loaded){
			store_motif_mapping = as_pot->get_motif_assignment();
			const double mm_choice = rand_gen->rand();
			if (mm_choice < motif_move_prob) {
				//as_pot->move_store_mapping();
				const bool as_result = as_pot->make_move(meta_data->get_pose());
				if (!as_result){
					std::cerr << "ca_motif_anneal_protocol: ERROR: as_geom_pot: move failure\n";
				}
			}
		}
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

	inline bool stage_rejected(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){
		state.num_rej++;
		state.tries++;
		//state.steps_no_update++;
		//TODO investigate more efficient alternative to this for big moves:
		state.was_prev_big_move = state.was_big_move;
		//protein->restore_backed_up_coords();
		protein->restore_coords(coords_backup);
		//as_pot->set_motif_assignment(store_motif_mapping);

		if (motif_loaded){
			//as_pot->move_restore_mapping();
			as_pot->set_motif_assignment(store_motif_mapping);
		}

		if (state.tries > this->recalc_tries){
			//meta_data->set_update_pair_lists_flag();
			state.steps_no_update = 0;
			state.cumulative_move = 0;
			meta_data->recalc_pair_lists_dists();
			state.energy = potentials->get_energy(meta_data, state.energy_map);
		    const double deltaE = state.energy - state.prev_acc_energy;
			if (fabs(deltaE) > energy_tolerance){
				state.prev_acc_energy = state.energy;
				std::cerr << "\nERROR: mc_protocol_interface: energy calculation problem\tdeltaE:\t" << deltaE << "\t"
						<< "at step:\t" << state.step_num
						<< std::endl;
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

	void stage_initialise_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){



		/*
			out_traj.open("traj.pdb", std::ios::out);
			std::cout << "bb_start\t" << state.step_num << "\t";
			state.energy_map.print_weighted_components(std::cout);

			protein->outputPdb(out_traj);
		 */


	}

	void stage_accepted_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){
		//std::cout << "bb_accepted\t" << state.step_num << "\t";
		//state.energy_map.print_weighted_components(std::cout);
		//protein->outputPdb(out_traj);
		//as_pot->print_assign_info(std::cout);
		//std::cout << std::endl;
	}

	void stage_rejected_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){
		//std::cout << "bb_rejected\t" << state.step_num << "\t";
		//state.energy_map.print_weighted_components(std::cout);
		//protein->outputPdb(out_traj);
		//as_pot->print_assign_info(std::cout);
		//std::cout << std::endl;
	}


};







}
}
}
}


#endif /* CA_MOTIF_ANNEAL_PROTOCOL_H_ */
