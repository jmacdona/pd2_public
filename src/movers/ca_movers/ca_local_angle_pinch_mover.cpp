/*
 * ca_local_angle_pinch_mover.cpp
 *
 *  Created on: 28 Sep 2010
 *      Author: jmacdona
 */

#include "ca_local_angle_pinch_mover.h"

using namespace PRODART::POSE;
using namespace PRODART::UTILS;

namespace PRODART {
namespace POSE {
namespace MOVERS {
namespace CA {


ca_local_angle_pinch_mover::ca_local_angle_pinch_mover(const double mdist,
		const int rnum){

	this->max_rot_ang = mdist;
	angleNum = rnum;

}

void ca_local_angle_pinch_mover::init(){

}

bool ca_local_angle_pinch_mover::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{
	mover_flags& return_flags = meta_data->get_mover_flags();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	const double angle_rot = this->rand_gen->rand(2.0*max_rot_ang) - max_rot_ang;
	//double move_dist = 0;


	double max_move = 0;


	const int originResNum = angleNum + 1;

	const vector3d origin =  protein->get_bb_atom_coords(POSE::CA, originResNum); // currentPdb->getVecByResidueIndex(originResNum, CA, chain);
	const vector3d origin_p1 = protein->get_bb_atom_coords(POSE::CA, originResNum+1); //currentPdb->getVecByResidueIndex(originResNum+1, CA, chain);
	const vector3d origin_m1 = protein->get_bb_atom_coords(POSE::CA, originResNum-1); //currentPdb->getVecByResidueIndex(originResNum-1, CA, chain);

	const vector3d rot_axis = (origin_p1 - origin) * (origin - origin_m1);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	const PRODART::POSE::const_chain_shared_ptr orig_chain = protein->get_residue(originResNum)->get_chain();


	{
		// minus 1
		const quaternion_rotate_params params(angle_rot, rot_axis_norm);
		atom_shared_ptr atm = protein->get_bb_atom(POSE::CA, originResNum-1);
		const vector3d startVec = atm->get_coords();
		vector3d endVec = startVec - origin;
		quaternion_rotate(endVec, params);
		endVec = endVec + origin;
		atm->set_coords(endVec);
		const double move_dist = (startVec - endVec).mod_sq();
		if (move_dist > max_move){
			max_move = move_dist;
		}
	}
	{
		// plus 1
		const quaternion_rotate_params params(-angle_rot, rot_axis_norm);
		atom_shared_ptr atm = protein->get_bb_atom(POSE::CA, originResNum+1);
		const vector3d startVec = atm->get_coords();
		vector3d endVec = startVec - origin;
		quaternion_rotate(endVec, params);
		endVec = endVec + origin;
		atm->set_coords(endVec);
		const double move_dist = (startVec - endVec).mod_sq();
		if (move_dist > max_move){
			max_move = move_dist;
		}
	}


	return_flags.move_completed = true;
	return_flags.move_dist = std::sqrt(max_move);

	return_flags.is_large_move = false;

	return true;

}
ca_local_angle_pinch_mover_shared_ptr new_ca_local_angle_pinch_mover(const double mdist,
			const int resNum){
	if (resNum >= 0){

		ca_local_angle_pinch_mover_shared_ptr ptr(new ca_local_angle_pinch_mover(mdist, resNum));
		return ptr;

	}
	std::cerr << "ERROR: new_ca_local_angle_pinch_mover: can not make mover\t" << resNum << std::endl;
	return ca_local_angle_pinch_mover_shared_ptr();
}

move_set_shared_ptr ca_local_angle_pinch_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_ang,
		const PRODART::POSE::MOVERS::bool_vector& allowed_residues){

	move_set_shared_ptr this_move_set = new_move_set();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();

	const int num_chains = protein->get_chain_count();

	for (int chain_num = 0; chain_num < num_chains; chain_num++){
		PRODART::POSE::chain_shared_ptr currChain = protein->get_chain(chain_num);
		const int start_resnum = currChain->get_first_internal_residue_index();
		const int end_resnum = currChain->get_last_internal_residue_index();
		for (int angNum = start_resnum; angNum <= end_resnum-2; angNum++){
			const int originResNum = angNum + 1;
			atom_shared_ptr o_m1_atm = protein->get_bb_atom(POSE::CA, originResNum-1);
			atom_shared_ptr o_atm = protein->get_bb_atom(POSE::CA, originResNum);
			atom_shared_ptr o_p1_atm = protein->get_bb_atom(POSE::CA, originResNum+1);

			if (o_m1_atm->isActive() && o_atm->isActive() && o_p1_atm->isActive()
					&& allowed_residues[originResNum-1] && allowed_residues[originResNum] && allowed_residues[originResNum+1]){
				mover_shared_ptr ptr = new_ca_local_angle_pinch_mover(max_ang,
						angNum);
				this_move_set->add_move(ptr, 1.0);

			}


		}
	}


	return this_move_set;
}

move_set_shared_ptr ca_local_angle_pinch_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_ang){
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	bool_vector allowed_residues(protein->get_residue_count(), true);
	return ca_local_angle_pinch_uni_dist_move_set_factory(meta_data,
			max_ang,
			allowed_residues);
}










}
}
}
}

