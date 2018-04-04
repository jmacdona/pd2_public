/*
 * structure_prediction.h
 *
 *  Created on: 9 Oct 2013
 *      Author: jmacdona
 */

#ifndef STRUCTURE_PREDICTION_H_
#define STRUCTURE_PREDICTION_H_

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "prodart_env/prodart_env.h"
#include <exception>
#include "pose_meta/ca_pose_meta.h"
#include "pose_meta/bb_pose_meta.h"
#include "pose_meta/frag_classify/fragment_classifier_interface.h"
#include "pose_meta/frag_classify/ca_fragment_classifier.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "utils/line_fit.h"
#include "pose/residue_type.h"
#include "pose/atom_type.h"
#include "pose/atom.h"
#include "pose/residue.h"
#include "pose/pose.h"
#include "pose_utils/pose_utils.h"
#include "prodart_env/prodart_env.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "movers/ca_movers/ca_single_dihedral_mover.h"
#include "movers/ca_movers/ca_single_angle_mover.h"
#include "movers/ca_movers/ca_single_crankshaft_mover.h"
#include "movers/ca_movers/ca_single_dihedral_mover.h"
#include "movers/ca_movers/ca_single_angle_mover.h"
#include "movers/ca_movers/ca_single_crankshaft_mover.h"
#include "movers/ca_movers/ca_local_bond_mover.h"
#include "movers/ca_movers/ca_local_angle_pinch_mover.h"
#include "movers/move_set_factory.h"
#include "pose_meta/ca_pose_meta.h"
#include "potentials/ca_potentials/ca_bump.h"
#include "potentials/bond_potential.h"
#include "potentials/potentials_container.h"
#include "potentials/potentials_factory.h"
#include "simulation/monte_carlo.h"
#include "simulation/minimiser.h"
#include "simulation/ca_mc_protocols/simple_ca_mc_protocol.h"
#include "simulation/ca_mc_protocols/simple_ca_sim_anneal_protocol.h"
#include "backbonebuilder/backbone_builder.h"
#include "potentials/potentials_container.h"
#include "simulation/bb_mc_protocols/fixed_bb_motif_anneal.h"
#include "simulation/mpi/mpi_replica_exchange_mc_protocol.h"
#include "simulation/mpi/mpi_replica_exchange_mc_protocol_verbose.h"


namespace PRODART {
namespace PROTOCOLS{
namespace MPI {
namespace STRUCTURE_PREDICTION {


//! predict structure using PSICOV prediction
bool mpi_replica_exchange_folding(PRODART::POSE::residue_type_vector seq,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const std::string psicov_file,
		const unsigned long REM_steps,
		boost::mpi::communicator comm,
		const std::string rep_specs,
		PRODART::POSE::pose_shared_ptr start_struct = PRODART::POSE::pose_shared_ptr(),
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");


}
}
}
}
#endif /* STRUCTURE_PREDICTION_H_ */
