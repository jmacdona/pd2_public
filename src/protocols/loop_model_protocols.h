/*
 * loop_model_protocols.h
 *
 *  Created on: 28 Sep 2010
 *      Author: jmacdona
 */

#ifndef LOOP_MODEL_PROTOCOLS_H_
#define LOOP_MODEL_PROTOCOLS_H_

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


namespace PRODART {
namespace PROTOCOLS{
namespace LOOP{

typedef std::vector<bool> bool_vector;
typedef std::map<POSE::atom_shared_ptr, UTILS::vector3d> atom_shared_ptr_vector3d_map;
typedef std::map<POSE::atom_shared_ptr, bool> atom_shared_ptr_bool_map;

PRODART::POSE::residue_shared_ptr_vector insert_blank_residues(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec);

PRODART::POSE::residue_shared_ptr_vector insert_ca_loop_linear_interpolate(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const double random_component = 0,
		const double CoM_component = 0.25);

PRODART::POSE::residue_shared_ptr_vector insert_ca_loop_linear_interpolate_minimise(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const double random_component = 0 );


PRODART::POSE::residue_shared_ptr_vector insert_ca_loop_linear_interpolate_minimise(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite);

PRODART::POSE::residue_shared_ptr_vector remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot);

PRODART::POSE::residue_shared_ptr_vector remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps = 500,
		const std::string ca_potentials_weight_set = "ca_default",
		const std::string bb_min_potentials_weight_set = "bb_min_default");

bool remodel_loops_multi(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long num_structs,
		std::ostream &output_structs,
		const bool renumber = false,
		const unsigned long min_steps = 100,
		const std::string ca_potentials_weight_set = "ca_default",
		const std::string bb_min_potentials_weight_set = "bb_min_default");

//! loop inserted before insert_resnum
bool remodel_loops_multi_filtered(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long num_structs,
		std::ostream &output_structs,
		const bool renumber = false,
		const unsigned long min_steps = 100,
		const std::string ca_potentials_weight_set = "ca_default",
		const std::string bb_min_potentials_weight_set = "bb_min_default",
		bool no_filter = false);

bool remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot);

bool remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps = 100,
		const std::string ca_potentials_weight_set = "ca_default",
		const std::string bb_min_potentials_weight_set = "bb_min_default");

bool run_ca_loop_anneal(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta);

bool run_ca_loop_anneal(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const std::string potentials_weight_set,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta);


bool ca_bb_multi_level_loop_design(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const bool_vector& as_site_mask,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta);

bool ca_bb_multi_level_loop_design(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta);

bool fixed_bb_motif_anneal(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_wt,
		const double final_wt);

bool fixed_bb_motif_brute_force_search(PRODART::POSE::pose_shared_ptr protein);

bool ca_minimise_bonds_bumps(PRODART::POSE::pose_shared_ptr protein,
		bool_vector& atom_selection);

bool ca_minimise_bonds_bumps(PRODART::POSE::pose_shared_ptr protein);

bool ca_minimise_bonds_bumps_angles(PRODART::POSE::pose_shared_ptr protein,
		bool_vector& atom_selection,
		bool auto_load_restraints = true);


double get_fiser_bb_loop_rmsd(PRODART::POSE::pose_shared_ptr ref,
		PRODART::POSE::pose_shared_ptr protein,
		const int start_res,
		const int length);

double get_fiser_bb_loop_rmsd_corrected(PRODART::POSE::pose_shared_ptr ref,
		PRODART::POSE::pose_shared_ptr protein,
		const int start_res,
		const int length);

bool_vector get_bb_hot_spots(PRODART::POSE::pose_shared_ptr protein);

void print_bb_hot_spots(PRODART::POSE::pose_shared_ptr protein, std::ostream& output);


inline atom_shared_ptr_vector3d_map get_coords(PRODART::POSE::pose_shared_ptr pose_,
		const bool_vector& atom_selection){
	atom_shared_ptr_vector3d_map rtn_vec;
	//rtn_vec.resize(atom_selection.size());
	for (unsigned int i = 0; i < atom_selection.size(); i++ ){
		if (atom_selection[i] == true){
			rtn_vec[pose_->get_atom(i)] = pose_->get_atom(i)->get_coords();
		}
	}
	return rtn_vec;
}

/*
inline atom_shared_ptr_vector3d_map get_backup_coords(PRODART::POSE::pose_shared_ptr pose_,
		const bool_vector& atom_selection){
	atom_shared_ptr_vector3d_map rtn_vec;
	//rtn_vec.resize(atom_selection.size());
	for (unsigned int i = 0; i < atom_selection.size(); i++ ){
		if (atom_selection[i] == true){
			rtn_vec[pose_->get_atom(i)] = pose_->get_atom(i)->get_backup_coords();
		}
	}

	return rtn_vec;
}
*/

/*
inline void restore_backup_coords(PRODART::POSE::pose_shared_ptr pose_,
		const atom_shared_ptr_vector3d_map& vec){

	atom_shared_ptr_vector3d_map::const_iterator it;

	for (it = vec.begin(); it != vec.end(); it++){
		it->first->set_backup_coords(it->second);
	}

}
*/

inline atom_shared_ptr_vector3d_map get_coords(const PRODART::POSE::atom_shared_ptr_vector& atom_selection){
	atom_shared_ptr_vector3d_map rtn_vec;
	//rtn_vec.resize(atom_selection.size());
	for (unsigned int i = 0; i < atom_selection.size(); i++ ){
		if (atom_selection[i]){
			rtn_vec[atom_selection[i]] = atom_selection[i]->get_coords();
		}
	}
	return rtn_vec;
}

inline void restore_coords(const atom_shared_ptr_vector3d_map& vec){

	atom_shared_ptr_vector3d_map::const_iterator it;

	for (it = vec.begin(); it != vec.end(); it++){
		it->first->set_coords(it->second);
	}

	/*
	for (unsigned int i = 0; i < atom_selection.size(); i++ ){
		if (atom_selection[i] == true){
			pose_->get_atom(i)->set_backup_coords(vec[i]);
		}
	}
	*/

}

atom_shared_ptr_bool_map get_activity_map(PRODART::POSE::pose_shared_ptr protein);

void restore_activity_map(const atom_shared_ptr_bool_map& vec);

void inactive_all_outside_radius(PRODART::POSE::pose_shared_ptr protein, const POSE::atom_shared_ptr_vector loop, const double radius );



}
}
}


#endif /* LOOP_MODEL_PROTOCOLS_H_ */
