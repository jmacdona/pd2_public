/*
 * mc_protocol_interface.h
 *
 *  Created on: 6 Aug 2010
 *      Author: jmacdona
 */

#ifndef MC_PROTOCOL_INTERFACE_H_
#define MC_PROTOCOL_INTERFACE_H_


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
#include <string>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "mc_protocol_functors.h"


namespace PRODART {
namespace POSE {
namespace SIM {

class mc_protocol_interface;

typedef boost::shared_ptr<mc_protocol_interface> mc_protocol_interface_shared_ptr;

//class mc_state;

class mc_protocol_interface {

private:


protected:



	//PRODART::POSE::POTENTIALS::potentials_energies_map energy_map;
	MTRand::MTRand_shared_ptr rand_gen;
	unsigned long steps_to_run;
	unsigned long max_tries;
	unsigned long recalc_tries;
	unsigned long max_steps_no_update;
	double beta;
	double energy_tolerance;
	double max_move_no_update;

	POSE::atom_shared_ptr_vector3d_map coords_backup;

	//std::ofstream traj_out;

	mc_protocol_interface(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr ptr)
	//: traj_out("temp_traj.pdb")
	{
		this->steps_to_run = steps;
		this->beta = _beta;
		max_tries = 10000;
		rand_gen = ptr;
		recalc_tries = 100;
		energy_tolerance = 0.00000000001;
		max_steps_no_update = 1000;
		max_move_no_update = 2;
	}


public:

	virtual ~mc_protocol_interface(){
	}

	void set_steps_to_run(const unsigned long  val){
		steps_to_run = val;
	}

	virtual bool stage_initialise(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	virtual void stage_initialise_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){}
	virtual bool stage_continue(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	virtual bool stage_move(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	virtual bool stage_accept_test(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	virtual bool stage_accepted(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	virtual void stage_accepted_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){}
	virtual bool stage_rejected(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	virtual void stage_rejected_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){}
	virtual bool stage_before_next_round(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	virtual bool stage_finish(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);




};





/*
 * INLINED functions
 */

inline bool mc_protocol_interface::stage_initialise(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){

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

	return true;

}
inline bool mc_protocol_interface::stage_continue(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
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
inline bool mc_protocol_interface::stage_move(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
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
inline bool mc_protocol_interface::stage_accept_test(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
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
inline bool mc_protocol_interface::stage_accepted(default_mc_state& state,
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

	/*
	std::cout << "accepted\t" << state.step_num << "\t";
	state.energy_map.print_weighted_components(std::cout);
	*/

	//protein->outputPdb(traj_out);

	return true;
}
inline bool mc_protocol_interface::stage_rejected(default_mc_state& state,
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
inline bool mc_protocol_interface::stage_before_next_round(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
	return true;
}
inline bool mc_protocol_interface::stage_finish(default_mc_state& state,
		const PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
		PRODART::POSE::MOVERS::mover_shared_ptr movers){
	state.final_pose = protein->clone();
	return true;
}







}
}
}


#endif /* MC_PROTOCOL_INTERFACE_H_ */
