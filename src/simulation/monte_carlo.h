/*
 * monte_carlo.h
 *
 *  Created on: 6 Aug 2010
 *      Author: jmacdona
 */

#ifndef MONTE_CARLO_H_
#define MONTE_CARLO_H_

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
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "mc_protocol_interface.h"
#include "ca_mc_protocols/simple_ca_mc_protocol.h"
#include "ca_mc_protocols/simple_ca_sim_anneal_protocol.h"
#include "ca_mc_protocols/ca_motif_anneal_protocol.h"
#include "bb_mc_protocols/simple_bb_mc_protocol.h"
#include "bb_mc_protocols/bb_motif_anneal_protocol.h"

namespace PRODART {
namespace POSE {
namespace SIM {


class monte_carlo;

typedef boost::shared_ptr<monte_carlo> monte_carlo_shared_ptr;

monte_carlo_shared_ptr new_monte_carlo(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use,
		PRODART::POSE::MOVERS::mover_shared_ptr movers_to_use,
		mc_protocol_interface_shared_ptr mc_protocol_to_use);

class monte_carlo : public PRODART::POSE::MOVERS::mover_interface {

	friend monte_carlo_shared_ptr new_monte_carlo(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use,
			PRODART::POSE::MOVERS::mover_shared_ptr movers_to_use,
			mc_protocol_interface_shared_ptr mc_protocol_to_use);

private:

	mutable PRODART::POSE::META::pose_meta_shared_ptr _pose_meta;
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr _potentials;
	PRODART::POSE::MOVERS::mover_shared_ptr _movers;
	mc_protocol_interface_shared_ptr _mc_protocol;



protected:

	//monte_carlo();
	monte_carlo(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use,
			PRODART::POSE::MOVERS::mover_shared_ptr movers_to_use,
			mc_protocol_interface_shared_ptr mc_protocol_to_use);

	mutable default_mc_state final_state;

public:


	bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const;


	default_mc_state get_state(){
		return this->final_state;
	}

};




}
}
}

#endif /* MONTE_CARLO_H_ */
