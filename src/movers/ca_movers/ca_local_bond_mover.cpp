/*
 * ca_local_bond_mover.cpp
 *
 *  Created on: 28 Sep 2010
 *      Author: jmacdona
 */
#include "ca_local_bond_mover.h"

using namespace PRODART::POSE;
using namespace PRODART::UTILS;

namespace PRODART {
namespace POSE {
namespace MOVERS {
namespace CA {



ca_local_bond_mover::ca_local_bond_mover(const double mdist,
		const int rnum){

	this->max_dist = mdist;
	bondNum = rnum;

}

void ca_local_bond_mover::init(){

}

bool ca_local_bond_mover::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{
	mover_flags& return_flags = meta_data->get_mover_flags();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	const double dist = this->rand_gen->rand(2.0*max_dist) - max_dist;
	//double move_dist = 0;

	const int orig_resnum =  bondNum;
	const int orig_resnum_p1 =  bondNum + 1;

	atom_shared_ptr o_atm = protein->get_bb_atom(POSE::CA, orig_resnum);
	atom_shared_ptr o_p1_atm = protein->get_bb_atom(POSE::CA, orig_resnum_p1);

	//if (o_atm && o_p1_atm){

	const vector3d o_coords = o_atm->get_coords();
	const vector3d o_p1_coords = o_p1_atm->get_coords();

	const vector3d move_vec = (o_p1_coords - o_coords) / (o_p1_coords - o_coords).mod();

	const vector3d new_o_coords = o_coords - (dist * move_vec);
	const vector3d new_o_p1_coords = o_p1_coords + (dist * move_vec);


	o_atm->set_coords(new_o_coords);
	o_p1_atm->set_coords(new_o_p1_coords);


	//}

	return_flags.move_completed = true;
	return_flags.move_dist = dist;

	return_flags.is_large_move = false;

	return true;

}
ca_local_bond_mover_shared_ptr new_ca_local_bond_mover(const double mdist,
			const int resNum){
	if (resNum >= 0){

		ca_local_bond_mover_shared_ptr ptr(new ca_local_bond_mover(mdist, resNum));
		return ptr;

	}
	std::cerr << "ERROR: new_ca_local_bond_mover: can not make mover\t" << resNum << std::endl;
	return ca_local_bond_mover_shared_ptr();
}

move_set_shared_ptr ca_local_bond_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_dist,
		const PRODART::POSE::MOVERS::bool_vector& allowed_residues){

	move_set_shared_ptr this_move_set = new_move_set();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();

	const int num_chains = protein->get_chain_count();

	for (int chain_num = 0; chain_num < num_chains; chain_num++){
		PRODART::POSE::chain_shared_ptr currChain = protein->get_chain(chain_num);
		const int start_resnum = currChain->get_first_internal_residue_index();
		const int end_resnum = currChain->get_last_internal_residue_index();
		for (int bondNum = start_resnum; bondNum <= end_resnum-1; bondNum++){
			const int orig_resnum =  bondNum;
			const int orig_resnum_p1 =  bondNum + 1;
			atom_shared_ptr o_atm = protein->get_bb_atom(POSE::CA, orig_resnum);
			atom_shared_ptr o_p1_atm = protein->get_bb_atom(POSE::CA, orig_resnum_p1);

			if (o_atm->isActive() && o_p1_atm->isActive()
					&& allowed_residues[orig_resnum] && allowed_residues[orig_resnum_p1]){
				mover_shared_ptr ptr = new_ca_local_bond_mover(max_dist,
						bondNum);
				this_move_set->add_move(ptr, 1.0);

			}


		}
	}


	return this_move_set;
}

move_set_shared_ptr ca_local_bond_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_dist){
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	bool_vector allowed_residues(protein->get_residue_count(), true);
	return ca_local_bond_uni_dist_move_set_factory(meta_data,
			max_dist,
			allowed_residues);
}



}
}
}
}


