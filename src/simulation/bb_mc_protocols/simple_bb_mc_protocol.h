/*
 * simple_bb_mc_protocol.h
 *
 *  Created on: 8 Nov 2010
 *      Author: jmacdona
 */

#ifndef SIMPLE_BB_MC_PROTOCOL_H_
#define SIMPLE_BB_MC_PROTOCOL_H_

#include "simulation/mc_protocol_interface.h"
#include <limits>
#include <iostream>
#include <fstream>

namespace PRODART {
namespace POSE {
namespace SIM {
namespace BB {


class simple_bb_mc_protocol;


typedef boost::shared_ptr<simple_bb_mc_protocol> simple_bb_mc_protocol_shared_ptr;

mc_protocol_interface_shared_ptr new_simple_bb_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr);

class simple_bb_mc_protocol : public mc_protocol_interface {


	friend mc_protocol_interface_shared_ptr new_simple_bb_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr);

private:

	std::ofstream out_traj;

protected:


	simple_bb_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr) : mc_protocol_interface(steps, _beta, mt_ptr){

	}

public:

	void stage_initialise_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){


		out_traj.open("traj.pdb", std::ios::out);

		std::cout << "bb_start\t" << state.step_num << "\t";
		state.energy_map.print_weighted_components(std::cout);

		protein->outputPdb(out_traj);

	}

	void stage_accepted_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){
		std::cout << "bb_accepted\t" << state.step_num << "\t";
		state.energy_map.print_weighted_components(std::cout);
		protein->outputPdb(out_traj);
	}

	void stage_rejected_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){
		std::cout << "bb_rejected\t" << state.step_num << "\t";
		state.energy_map.print_weighted_components(std::cout);
		//protein->outputPdb(out_traj);
	}


};







}
}
}
}

#endif /* SIMPLE_BB_MC_PROTOCOL_H_ */
