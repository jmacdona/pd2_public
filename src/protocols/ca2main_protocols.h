/*
 * ca2main_protocols.h
 *
 *  Created on: 24 Sep 2010
 *      Author: jmacdona
 */

#ifndef CA2MAIN_PROTOCOLS_H_
#define CA2MAIN_PROTOCOLS_H_

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
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
#include "pose_utils/sse_axes.h"
#include "prodart_env/prodart_env.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "movers/ca_movers/ca_single_dihedral_mover.h"
#include "movers/ca_movers/ca_single_angle_mover.h"
#include "movers/ca_movers/ca_single_crankshaft_mover.h"
#include "pose_meta/ca_pose_meta.h"
#include "potentials/ca_potentials/ca_bump.h"
#include "potentials/bond_potential.h"
#include "potentials/potentials_container.h"
#include "potentials/potentials_factory.h"
#include "simulation/monte_carlo.h"
#include "simulation/ca_mc_protocols/simple_ca_mc_protocol.h"
#include "simulation/minimiser.h"
#include "backbonebuilder/backbone_builder.h"
#include "backbonebuilder/alphabet_bb_builder.h"
#include "backbonebuilder/alphabet_bb_builder_manager.h"

namespace PRODART {
namespace PROTOCOLS{
namespace CA2MAIN{

typedef std::vector<bool> bool_vector;
typedef std::vector<int> int_vector;
typedef std::vector<int_vector> int_vector_vector;
typedef std::map<POSE::atom_shared_ptr, POSE::atom_shared_ptr> atom_pair_map;
typedef std::map<POSE::residue_shared_ptr, POSE::residue_shared_ptr> residue_pair_map;


int_vector_vector get_strands_vector(const PRODART::POSE::three_state_sec_struct_vector& secs);
int_vector_vector get_helices_vector(const PRODART::POSE::three_state_sec_struct_vector& secs);
int_vector_vector get_loops_vector(const PRODART::POSE::three_state_sec_struct_vector& secs);

//! add 2 end atoms (in new residues) to the terminii of each chain (if isPeptide)
bool add_two_end_ca_atoms_residues(PRODART::POSE::pose_shared_ptr protein);
//! delete 2 end atoms (and containing residues) from the terminii of each chain (if isPeptide)
bool delete_two_end_residues(PRODART::POSE::pose_shared_ptr protein);

// simple ca2main with no CA refinement, minimisation or optimisation
bool simple_ca2main(PRODART::POSE::pose_shared_ptr protein);

// simple ca2main with no CA refinement, minimisation or optimisation
bool simple_ca2main(PRODART::POSE::pose_shared_ptr protein, const bool_vector& loop_mask);

// alphabet ca2main with no CA refinement, minimisation or optimisation
bool new_ca2main(PRODART::POSE::pose_shared_ptr protein);

//! alphabet ca2main with no CA refinement, minimisation or optimisation - NOTE that loop_mask is assumed to be the moved CA atoms so that peptide bonds will be rebuilt between start-1 and end+1
bool new_ca2main(PRODART::POSE::pose_shared_ptr protein, const bool_vector& loop_mask);

// simple ca2main with minimisation
bool simple_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

// alphabet ca2main with minimisation
bool new_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

// alphabet ca2main with minimisation
bool new_ca2main_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

// alphabet ca2main with minimisation
bool new_ca2main_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

// alphabet ca2main with minimisation
bool simple_ca2main_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

// alphabet ca2main with minimisation
bool new_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

// alphabet ca2main with minimisation
bool new_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot);

bool ca_mc_refine(PRODART::POSE::pose_shared_ptr protein, const unsigned long steps, const double beta, const std::string ca_pot);

// simple ca2main with sec struct refinement and minimisation
bool simple_ca2main_secs_refine_minimise(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const unsigned long ca_ref_steps,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

// simple ca2main with minimisation
bool simple_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot);

// simple ca2main with minimisation
bool simple_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");

bool bb_minimise(PRODART::POSE::pose_shared_ptr protein,
		const std::string pot_preset,
		const unsigned long min_steps = 100);
bool bb_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot,
		const unsigned long min_steps = 100);
bool bb_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const std::string pot_preset,
		const unsigned long min_steps = 100);

bool bb_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const std::string pot_preset,
		const unsigned long min_steps = 100);
bool bb_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot,
		const unsigned long min_steps = 100);
bool bb_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const std::string pot_preset,
		const unsigned long min_steps = 100);

bool bb_minimise_geom(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps = 100);
bool bb_minimise_geom(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps = 100);

PRODART::POSE::atom_shared_ptr_vector get_loop_moved_atoms(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& residue_loop_mask);
PRODART::POSE::atom_shared_ptr_vector get_loop_moved_atoms(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum);


PRODART::POSE::atom_shared_ptr_vector get_loop_moved_N_CA_C_O_atoms(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& residue_loop_mask);
PRODART::POSE::atom_shared_ptr_vector get_loop_moved_N_CA_C_O_atoms(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum);



bool_vector make_residue_loop_mask(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum);

bool_vector make_residue_loop_mask(const PRODART::POSE::three_state_sec_struct_vector& secs);

// simple ca2main with sec struct refinement and minimisation of loops and sec struct independently
bool simple_ca2main_secs_indep_refine_minimise(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");


// ca2main with no CA refinement, minimisation or optimisation with new alphabet method
bool alphabet_ca2main_fixed_ca_minimise(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::BB_BUILDER::alphabet_bb_builder_shared_ptr bbb,
		const unsigned long min_steps = 500,
		const std::string bb_min_potentials_weight_set = "bb_min_default");


//!get RMSD for ca2main method - special for Ben's ca2main paper. Ignore first NCONCON and last CONCONCO
double ca2main_get_rmsd_special(PRODART::POSE::pose_shared_ptr ref, PRODART::POSE::pose_shared_ptr protein);

//!get RMSD for ca2main method - special for Ben's ca2main paper. sidechain only (not including CB)
double ca2main_get_rmsd_sidechain(PRODART::POSE::pose_shared_ptr ref, PRODART::POSE::pose_shared_ptr protein);



bool auto_add_secs_restraints(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const double weight);

}
}
}


#endif /* CA2MAIN_PROTOCOLS_H_ */
