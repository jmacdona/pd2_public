//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * ca_single_dihedral_mover.cpp
 *
 *  Created on: 2 Mar 2010
 *      Author: jmacdona
 */

#include "ca_single_dihedral_mover.h"



namespace PRODART {
namespace POSE {
namespace MOVERS {
namespace CA {

using namespace PRODART::POSE;
using namespace PRODART::UTILS;

ca_single_dihedral_uni_dist_mover_shared_ptr new_ca_single_dihedral_uni_dist_mover(const bool is_forwards,
		const double max_angle,
		const int ang_num,
		const int largest_move_resnum){

	ca_single_dihedral_uni_dist_mover_shared_ptr ptr(new ca_single_dihedral_uni_dist_mover(is_forwards, max_angle, ang_num, largest_move_resnum));
	return ptr;

}

move_set_shared_ptr ca_single_dihedral_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data, const double max_angle){
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	bool_vector allowed_residues(protein->get_residue_count(), true);
	return ca_single_dihedral_uni_dist_move_set_factory(meta_data, max_angle, allowed_residues);
}

move_set_shared_ptr ca_single_dihedral_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data, const double max_angle, const bool_vector& allowed_residues){
	move_set_shared_ptr this_move_set = new_move_set();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	//protein->index();

	const int num_chains = protein->get_chain_count();

	for (int chain_num = 0; chain_num < num_chains; chain_num++){
		PRODART::POSE::chain_shared_ptr currChain = protein->get_chain(chain_num);
		const int start_resnum = currChain->get_first_internal_residue_index();
		const int end_resnum = currChain->get_last_internal_residue_index();
		const int half_way = ((end_resnum - start_resnum + 1) / 2) + start_resnum - 3;
		for (int ang_num = start_resnum; ang_num <= end_resnum - 3; ang_num++){
			if (ang_num > half_way){
				// forwards
				bool allowed = true;
				const int originResNum = ang_num + 1;
				int i = 0 , bm_resnum = 0;
				for ( i = originResNum+1; i <= end_resnum; i++){
					atom_shared_ptr atm = protein->get_bb_atom(PRODART::POSE::CA,i);
					if (atm->isActive() == false || allowed_residues[i] == false) {
						if (allowed_residues[i] == false) allowed = false;
						break;
					}
				}
				bm_resnum = i-1;
				//std::cout << originResNum << "\t" << bm_resnum << "\t" << end_resnum << "\t" << std::endl;
				if (bm_resnum != originResNum && allowed){
					mover_shared_ptr ptr = new_ca_single_dihedral_uni_dist_mover(true,
							max_angle,
							ang_num,
							bm_resnum);
					this_move_set->add_move(ptr, 1.0);
				}
			}
			else {
				// backwards
				bool allowed = true;
				const int originResNum = ang_num + 2;
				int i = 0 , bm_resnum = 0;
				for ( i = originResNum -1; i >= start_resnum; i--){
					atom_shared_ptr atm = protein->get_bb_atom(PRODART::POSE::CA,i);
					if (atm->isActive() == false || allowed_residues[i] == false) {
						//bm_resnum = 1;
						if (allowed_residues[i] == false) allowed = false;
						break;
					}
				}
				bm_resnum = i+1;
				//std::cout << originResNum << "\t" << bm_resnum << "\t" << start_resnum << "\t" << std::endl;
				if (bm_resnum != originResNum && allowed){
					mover_shared_ptr ptr = new_ca_single_dihedral_uni_dist_mover(false,
							max_angle,
							ang_num,
							bm_resnum);
					this_move_set->add_move(ptr, 1.0);
				}
			}
		}

	}

	return this_move_set;
}

bool ca_single_dihedral_uni_dist_mover::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const {
		//mover_flags& return_flags) const{

	mover_flags& return_flags = meta_data->get_mover_flags();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();

	//const vector3d orig_coords = protein->get_bb_atom_coords(PRODART::POSE::CA,this->big_move_resnum);

	const double ang_rot = this->rand_gen->rand(2.0*max_rot_ang) - max_rot_ang;
	double move_dist = 0;
	if (isFwds){
		move_dist = PRODART::POSE_UTILS::KINE::rotate_dihedral_forwards_ca_track_move(protein, angle_num, ang_rot);
	}
	else {
		move_dist = PRODART::POSE_UTILS::KINE::rotate_dihedral_backwards_ca_track_move(protein, angle_num, ang_rot);
	}

	//const vector3d end_coords = protein->get_bb_atom_coords(PRODART::POSE::CA,this->big_move_resnum);


	return_flags.move_completed = true;
	//TODO sort out how to check for large moves
	return_flags.move_dist = move_dist;//(orig_coords - end_coords).mod();

/*
	std::cout << angle_num << "\t"
			<< this->big_move_resnum << "\t"
			<< ang_rot << "\t"
			<< return_flags.move_dist << std::endl;
*/

	return_flags.is_large_move = false;

	return true;
}



ca_single_dihedral_uni_dist_mover::ca_single_dihedral_uni_dist_mover() : isFwds(true),
		max_rot_ang(0), angle_num(0), big_move_resnum(0) {
	this->init();
}

ca_single_dihedral_uni_dist_mover::ca_single_dihedral_uni_dist_mover(const bool is_forwards,
		const double max_angle,
		const int ang_num,
		const int largest_move_resnum){
	this->init();
	isFwds = is_forwards;
	max_rot_ang = max_angle;
	angle_num = ang_num;
	this->big_move_resnum = largest_move_resnum;
}

void ca_single_dihedral_uni_dist_mover::init(){
	isFwds = true;
	max_rot_ang = 0;
	angle_num = 0;
	this->type = mover_interface::mt_pure_mover;
}



}
}
}
}

