//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose_meta_interface.cpp
 *
 *  Created on: 5 Mar 2010
 *      Author: jmacdona
 */

#include "pose_meta_interface.h"



using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;

namespace PRODART {
namespace POSE {
namespace META {




pose_meta_interface::pose_meta_interface(PRODART::POSE::pose_shared_ptr pose_to_attach){

	this->pose_ = pose_to_attach;
	pose_->index();
	refresh_pair_lists = true;
	this->gradient.clear();
	this->gradient.resize(pose_->get_all_atom_count(), PRODART::UTILS::vector3d(0,0,0));
	this->residue_conformation_class.resize(pose_->get_residue_count(), PRODART::POSE::ss4_OTHER);
	this->residue_sec_struct.resize(pose_->get_residue_count(), PRODART::POSE::ss3_OTHER);
	sec_struct_restraint_list.resize(pose_->get_residue_count(), PRODART::POSE::ss4_UNDEF);
	sec_struct_restraint_weight_list.resize(pose_->get_residue_count(), 0.0);

	strand_hbond_count.clear();
	strand_hbond_count.resize(pose_->get_residue_count(), 0);
	helix_hbond_count.clear();
	helix_hbond_count.resize(pose_->get_residue_count(), 0);

	phi_psi_omega_vec.clear();
	phi_psi_omega_vec.resize(pose_->get_residue_count(), double_3_tuple(PI, PI, PI));

	move_margin = 2.0;

	this->clear_ca_GO_restraints();
	go_restraints_set = false;
	//this->updateBondedPairList(this->bond_list);




}

void pose_meta_interface::make_update_state_store(){
	update_state_store.clear();
	const int num_atms = pose_->get_all_atom_count();

	for (int i = 0; i < num_atms; i++){
		atom_shared_ptr atm = pose_->get_atom(i);
		update_state_store[atm] = atm->get_coords();
	}

}

void pose_meta_interface::make_has_moved_coord_store(){
	has_moved_coord_store.clear();
	const int num_atms = pose_->get_all_atom_count();

	for (int i = 0; i < num_atms; i++){
		atom_shared_ptr atm = pose_->get_atom(i);
		has_moved_coord_store[atm] = atm->get_coords();
	}
}

bool  pose_meta_interface::pair_lists_ok(){
	std::map<POSE::atom_shared_ptr, PRODART::UTILS::vector3d>::const_iterator it;
	const double max_distsq_change = move_margin * move_margin;
	for (it = update_state_store.begin(); it != update_state_store.end(); it++ ){
		double distsq_change = (it->first->get_coords() - it->second).mod_sq();
		if (distsq_change > max_distsq_change){
			return false;
		}
	}
	return true;
}

void pose_meta_interface::reset_update_state_store(){
	std::map<POSE::atom_shared_ptr, PRODART::UTILS::vector3d>::iterator it;
	for (it = update_state_store.begin(); it != update_state_store.end(); it++ ){
		it->second = it->first->get_coords();
	}
}

void pose_meta_interface::reset_has_moved_coord_store(){
	std::map<POSE::atom_shared_ptr, PRODART::UTILS::vector3d>::iterator it;
	for (it = has_moved_coord_store.begin(); it != has_moved_coord_store.end(); it++ ){
		it->second = it->first->get_coords();
	}
}

void pose_meta_interface::update_has_moved_flags(){
	std::map<POSE::atom_shared_ptr, PRODART::UTILS::vector3d>::iterator it;
	for (it = has_moved_coord_store.begin(); it != has_moved_coord_store.end(); it++ ){
		const vector3d stored_coords = it->second;
		const vector3d actual_coords = it->first->get_coords();
		if (stored_coords == actual_coords){
			it->first->set_has_moved_flag(false);
		}
		else {
			it->first->set_has_moved_flag(true);
		}
	}
}

void pose_meta_interface::updateBondedPairList(bonded_ele_vector& in_bond_list){
	bonded_ele_vector::iterator it;

	for (it = in_bond_list.begin(); it != in_bond_list.end(); it++){

		if (it->atom1_ptr
				&& it->atom2_ptr){
			if (it->atom1_ptr->isActive()
					&& it->atom2_ptr->isActive()){
				const vector3d v1 = it->atom1_ptr->get_coords();
				const vector3d v2 = it->atom2_ptr->get_coords();

				const double dist_sq = (v1 - v2).mod_sq();
				const double dist =std::sqrt(dist_sq);

				it->dist_sq = dist_sq;
				it->dist = dist;
			}
		}
		else {
			//cerr << "WTF!\t" << it->atom1 << endl;
		}

	}
}

void pose_meta_interface::updateCoordRestraintList(coord_rst_element_vector& coord_rst_list){
	coord_rst_element_vector::iterator it;
	//cerr << "updateCoordRestraintList: updating: " << coord_rst_list.size() << endl;
	for (it = coord_rst_list.begin(); it != coord_rst_list.end(); it++){

		if (it->atom_ptr){
			if (it->atom_ptr->isActive()){
				const vector3d v1 = it->atom_ptr->get_coords();
				const vector3d v2 = it->coord;

				const double dist_sq = (v1 - v2).mod_sq();
				const double dist =std::sqrt(dist_sq);

				it->dist_sq = dist_sq;
				it->dist = dist;
			}
			else {
				//cerr << "WTF!\t" << it->atom_ptr->get_type().get_label() << endl;
			}
		}
		else {
			//cerr << "WTF!\t" << it->atom_ptr->get_type().get_label() << endl;
			//cerr << "WTF!\t" << it->atom << endl;
		}

	}
}

void pose_meta_interface::updateNonBondedPairList(nb_ele_vector& bond_list){
	nb_ele_vector::iterator it;

	for (it = bond_list.begin(); it != bond_list.end(); it++){

		if (it->atom1_ptr
				&& it->atom2_ptr){
			if (it->atom1_ptr->isActive()
					&& it->atom2_ptr->isActive()){
				const vector3d v1 = it->atom1_ptr->get_coords();
				const vector3d v2 = it->atom2_ptr->get_coords();

				const double dist_sq = (v1 - v2).mod_sq();
				const double dist =std::sqrt(dist_sq);

				it->dist_sq = dist_sq;
				it->dist = dist;
			}
		}
		else {
			//cerr << "WTF!\t" << it->atom1 << endl;
		}

	}

}

void pose_meta_interface::updateAngleList(angle_ele_vector& angle_list){

	angle_ele_vector::iterator it;

	for (it = angle_list.begin(); it != angle_list.end(); it++){

		if (it->atom1_ptr
				&& it->atom2_ptr
				&& it->atom3_ptr){
			if (it->atom1_ptr->isActive()
					&& it->atom2_ptr->isActive()
					&& it->atom3_ptr->isActive()){
				const vector3d v1 = it->atom1_ptr->get_coords();
				const vector3d v2 = it->atom2_ptr->get_coords();
				const vector3d v3 = it->atom3_ptr->get_coords();

				const double ang = angle(v1, v2, v3);

				it->angle = ang;
			}
		}
		else {
			cerr << "WTF!\t" << it->atom1 << endl;
		}

	}


}

void pose_meta_interface::updateDihedralList(dihedral_ele_vector& dih_list){
	dihedral_ele_vector::iterator it;

	for (it = dih_list.begin(); it != dih_list.end(); it++){

		if (it->atom1_ptr
				&& it->atom2_ptr
				&& it->atom3_ptr
				&& it->atom4_ptr){
			if (it->atom1_ptr->isActive()
					&& it->atom2_ptr->isActive()
					&& it->atom3_ptr->isActive()
					&& it->atom4_ptr->isActive()){
				const vector3d v1 = it->atom1_ptr->get_coords();
				const vector3d v2 = it->atom2_ptr->get_coords();
				const vector3d v3 = it->atom3_ptr->get_coords();
				const vector3d v4 = it->atom4_ptr->get_coords();

				const double ang = dihedral(v1, v2, v3, v4);

				it->dih_angle = ang;
			}
		}
		else {
			cerr << "WTF!\t" << it->atom1 << endl;
		}

	}
}

void pose_meta_interface::updateDihedralList(simple_harmonic_dihedral_element_vector& dih_list){
	simple_harmonic_dihedral_element_vector::iterator it;

	for (it = dih_list.begin(); it != dih_list.end(); it++){

		if (it->atom1_ptr
				&& it->atom2_ptr
				&& it->atom3_ptr
				&& it->atom4_ptr){
			if (it->atom1_ptr->isActive()
					&& it->atom2_ptr->isActive()
					&& it->atom3_ptr->isActive()
					&& it->atom4_ptr->isActive()){
				const vector3d v1 = it->atom1_ptr->get_coords();
				const vector3d v2 = it->atom2_ptr->get_coords();
				const vector3d v3 = it->atom3_ptr->get_coords();
				const vector3d v4 = it->atom4_ptr->get_coords();

				const double ang = dihedral(v1, v2, v3, v4);

				it->dih_angle = ang;
			}
		}
		else {
			cerr << "WTF!\t" << it->atom1 << endl;
		}

	}
}

void pose_meta_interface::add_sse_restraint(const int start_resnum,
		const int end_resnum,
		POSE::four_state_sec_struct sse_type,
		const double half_bond_const,
		const UTILS::vector3d start,
		const UTILS::vector3d end){

	const int atom_count = pose_->get_all_atom_count();
	const int res_count = pose_->get_residue_count();


	if (start_resnum < end_resnum){
		if (start_resnum < res_count && end_resnum < res_count){

			POSE::residue_shared_ptr res1 = this->pose_->get_residue(start_resnum);
			POSE::residue_shared_ptr res2 = this->pose_->get_residue(end_resnum);

			if (res1 && res2){

				POSE::atom_shared_ptr at1 = res1->get_bb_atom(POSE::CA);
				POSE::atom_shared_ptr at2 = res2->get_bb_atom(POSE::CA);

				if (at1 && at2){

					const int at1_seq_num = at1->get_seq_num();
					const int at2_seq_num = at2->get_seq_num();

					if (at1_seq_num < atom_count && at2_seq_num < atom_count){
						sse_axis_element ele;

						ele.start.coord = start;
						ele.end.coord = end;

						ele.start.res_num = start_resnum;
						ele.end.res_num = end_resnum;
						ele.start.atom_ptr = at1;
						ele.end.atom_ptr = at2;
						ele.start.atom = at1->get_seq_num();
						ele.end.atom = at2->get_seq_num();
						ele.secs = sse_type;
						ele.seg_intersect_dist_half_bond_const = half_bond_const;
						ele.seg_intersect_dih_half_bond_const = half_bond_const * 10;
						ele.seg_intersect_ang_half_bond_const = half_bond_const * 10;
						//ele.equilibrium_dist = eq_dist;
						//ele.half_bond_const = half_bond_const;
						//ele.dist_sq = (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()).mod_sq();
						//ele.dist = std::sqrt(ele.dist_sq);
						this->sse_axis_restraint_list.push_back(ele);


						if ((!at1->isActive()) || (!at2->isActive())){
							std::cerr << "WARNING: pose_meta_interface: add_sse_restraint: atoms not active so will have no effect" << std::endl;
						}

					}
					else{
						std::cerr << "ERROR: pose_meta_interface: can't add_sse_restraint: sequence numbers not valid" << std::endl;
					}
				}
				else{
					std::cerr << "ERROR: pose_meta_interface: can't add_sse_restraint: atom pointers are NULL" << std::endl;
				}
			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add_sse_restraint: residue pointers are not valid" << std::endl;
			}
		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add_sse_restraint: residue numbers not valid" << std::endl;
		}
	}
	else{
		std::cerr << "ERROR: pose_meta_interface: can't add_sse_restraint: residue numbers not valid - start_resnum must be less than end_resnum" << std::endl;
	}
}




}
}
}
