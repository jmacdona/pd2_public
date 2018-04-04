/*
 * ca_single_crankshaft_mover.cpp
 *
 *  Created on: 8 Sep 2010
 *      Author: jmacdona
 */

#include "ca_single_crankshaft_mover.h"




namespace PRODART {
namespace POSE {
namespace MOVERS {
namespace CA {

using namespace PRODART::POSE;
using namespace PRODART::UTILS;

ca_single_crankshaft_uni_dist_mover_shared_ptr new_ca_single_crankshaft_uni_dist_mover(const double max_angle,
			const int resNum1,
			const int resNum2){

	if (resNum2 > resNum1){

		ca_single_crankshaft_uni_dist_mover_shared_ptr ptr(new ca_single_crankshaft_uni_dist_mover(max_angle, resNum1, resNum2));
		return ptr;

	}
	std::cerr << "ERROR: new_ca_single_crankshaft_uni_dist_mover: can not make mover\t" << resNum1 << "\t" << resNum2 << std::endl;
	return ca_single_crankshaft_uni_dist_mover_shared_ptr();
}


ca_single_crankshaft_uni_dist_mover::ca_single_crankshaft_uni_dist_mover(const double max_angle,
		const int r1,
		const int r2){
	this->init();
	max_rot_ang = max_angle;
	resNum1 = r1;
	resNum2 = r2;
}

void ca_single_crankshaft_uni_dist_mover::init(){
	max_rot_ang = 0;
	resNum1 = 0;
	resNum2 = 0;
}


bool ca_single_crankshaft_uni_dist_mover::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{
	mover_flags& return_flags = meta_data->get_mover_flags();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();

	//const vector3d orig_coords = protein->get_bb_atom_coords(PRODART::POSE::CA,this->big_move_resnum);

	const double ang_rot = this->rand_gen->rand(2.0*max_rot_ang) - max_rot_ang;
	double move_dist = 0;
	move_dist = PRODART::POSE_UTILS::KINE::crankshaft_ca_track_move(protein, resNum1, resNum2, ang_rot);


	//const vector3d end_coords = protein->get_bb_atom_coords(PRODART::POSE::CA,this->big_move_resnum);


	return_flags.move_completed = true;
	return_flags.move_dist = move_dist;//(orig_coords - end_coords).mod();

/*
	std::cout << this->resNum1 << "\t"
			<< this->resNum2 << "\t"
			<< ang_rot << "\t"
			<< return_flags.move_dist << std::endl;
*/

	return_flags.is_large_move = false;

	return true;
}

move_set_shared_ptr ca_single_crankshaft_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_angle,
		const int min_seq_sep,
		const int max_seq_sep){
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	bool_vector allowed_residues(protein->get_residue_count(), true);
	return ca_single_crankshaft_uni_dist_move_set_factory(meta_data,
			max_angle,
			min_seq_sep,
			max_seq_sep,
			allowed_residues);
}


move_set_shared_ptr ca_single_crankshaft_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_angle,
		const int min_seq_sep,
		const int max_seq_sep,
		const PRODART::POSE::MOVERS::bool_vector& allowed_residues){
	move_set_shared_ptr this_move_set = new_move_set();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();

	const int num_chains = protein->get_chain_count();

	for (int chain_num = 0; chain_num < num_chains; chain_num++){
		PRODART::POSE::chain_shared_ptr currChain = protein->get_chain(chain_num);
		const int start_resnum = currChain->get_first_internal_residue_index();
		const int end_resnum = currChain->get_last_internal_residue_index();
		for (int resNum = start_resnum; resNum <= end_resnum; resNum++){
			for (int seq_sep = min_seq_sep ; seq_sep <= max_seq_sep; seq_sep++){

				bool allOK = true;
				for (int i = resNum+1; i < resNum + seq_sep; i++){
					if (i <= end_resnum){
						atom_shared_ptr atm = protein->get_bb_atom(PRODART::POSE::CA,i);
						if ( atm->isActive() == false || allowed_residues[i] == false) allOK = false;
					}
					else {
						allOK = false;
					}
				}

				if (resNum + seq_sep > end_resnum){
					allOK = false;
				}

				if (allOK){
					mover_shared_ptr ptr = new_ca_single_crankshaft_uni_dist_mover(max_angle,
							resNum,
							resNum+seq_sep);
					this_move_set->add_move(ptr, 1.0);
				}

			}


		}

	}


	return this_move_set;
}



}
}
}
}


