/*
 * monte_carlo.cpp
 *
 *  Created on: 6 Aug 2010
 *      Author: jmacdona
 */

#include "monte_carlo.h"

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

namespace PRODART {
namespace POSE {
namespace SIM {


monte_carlo_shared_ptr new_monte_carlo(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use,
		PRODART::POSE::MOVERS::mover_shared_ptr movers_to_use,
		mc_protocol_interface_shared_ptr mc_protocol_to_use){
	monte_carlo_shared_ptr ptr(new monte_carlo(potentials_to_use, movers_to_use, mc_protocol_to_use));
	return ptr;
}

bool monte_carlo::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{

	_pose_meta = meta_data;
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();

	default_mc_state state;

	_mc_protocol->stage_initialise(state,
			protein,
			_pose_meta,
			_potentials,
			_movers);

	_mc_protocol->stage_initialise_print(state,
			protein,
			_pose_meta,
			_potentials,
			_movers);

	while (_mc_protocol->stage_continue(state,
			protein,
			_pose_meta,
			_potentials,
			_movers) && !state.stop){

		_mc_protocol->stage_move(state,
					protein,
					_pose_meta,
					_potentials,
					_movers);

		if (_mc_protocol->stage_accept_test(state,
				protein,
				_pose_meta,
				_potentials,
				_movers)){
			_mc_protocol->stage_accepted(state,
							protein,
							_pose_meta,
							_potentials,
							_movers);
			_mc_protocol->stage_accepted_print(state,
							protein,
							_pose_meta,
							_potentials,
							_movers);
		}
		else {
			_mc_protocol->stage_rejected_print(state,
							protein,
							_pose_meta,
							_potentials,
							_movers);
			_mc_protocol->stage_rejected(state,
							protein,
							_pose_meta,
							_potentials,
							_movers);

		}

		_mc_protocol->stage_before_next_round(state,
						protein,
						_pose_meta,
						_potentials,
						_movers);

	}


	const bool final_result = _mc_protocol->stage_finish(state,
				protein,
				_pose_meta,
				_potentials,
				_movers);

	this->final_state = state;

	return final_result;
}


monte_carlo::monte_carlo( PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use,
		PRODART::POSE::MOVERS::mover_shared_ptr movers_to_use,
		mc_protocol_interface_shared_ptr mc_protocol_to_use) :
		_potentials(potentials_to_use), _movers(movers_to_use), _mc_protocol(mc_protocol_to_use){


}




}
}
}
