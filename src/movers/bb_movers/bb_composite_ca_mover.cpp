/*
 * bb_composite_ca_mover.cpp
 *
 *  Created on: 4 Nov 2010
 *      Author: jmacdona
 */
#include "bb_composite_ca_mover.h"
#include "bb_movers.h"

using namespace PRODART::POSE;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::PROTOCOLS;
using namespace PRODART::POSE::POTENTIALS;

using namespace std;

namespace PRODART {
namespace POSE {
namespace MOVERS {
namespace BB {




bool bb_composite_ca_mover::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{
	mover_flags& return_flags = meta_data->get_mover_flags();
	pose_shared_ptr pose_ = meta_data->get_pose();

	bool_vector loop_mask(pose_->get_residue_count(), false);

	residue_shared_ptr_vector::const_iterator iter;
	for (int i = this->resNum1; i <= this->resNum2; i++){
		loop_mask[i] = true;
	}

	bool_vector atom_selection(pose_->get_all_atom_count(), false);
	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){

			if (i > 0){
				if (pose_->get_bb_atom(C, i-1)->isActive()) atom_selection[pose_->get_bb_atom(C, i-1)->get_seq_num()] = true;
				if (pose_->get_bb_atom(O, i-1)->isActive()) atom_selection[pose_->get_bb_atom(O, i-1)->get_seq_num()] = true;
			}

			if (pose_->get_bb_atom(N, i)->isActive()) atom_selection[pose_->get_bb_atom(N, i)->get_seq_num()] = true;
			if (pose_->get_bb_atom(POSE::CA, i)->isActive()) atom_selection[pose_->get_bb_atom(POSE::CA, i)->get_seq_num()] = true;
			if (pose_->get_bb_atom(C, i)->isActive()) atom_selection[pose_->get_bb_atom(C, i)->get_seq_num()] = true;
			if (pose_->get_bb_atom(CB, i)->isActive()) atom_selection[pose_->get_bb_atom(CB, i)->get_seq_num()] = true;
			if (pose_->get_bb_atom(O, i)->isActive()) atom_selection[pose_->get_bb_atom(O, i)->get_seq_num()] = true;
			if (pose_->get_bb_atom(H, i)->isActive()) atom_selection[pose_->get_bb_atom(H, i)->get_seq_num()] = true;
			if (i+1 < loop_mask.size()){
				if (pose_->get_bb_atom(N, i+1)->isActive()) atom_selection[pose_->get_bb_atom(N, i+1)->get_seq_num()] = true;
				if (pose_->get_bb_atom(H, i+1)->isActive()) atom_selection[pose_->get_bb_atom(H, i+1)->get_seq_num()] = true;
			}
		}
	}

	LOOP::atom_shared_ptr_vector3d_map stored_backup = LOOP::get_coords(pose_, atom_selection);

	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container(this->ca_pot_set);
	if (ca_as_geom_pot_changed){
		//cout << "ca " << ca_as_geom_pot_wt << endl;
		ca_pot->get_default_energies_map().adjust_weight(potentials_name("as_geom_pot"), ca_as_geom_pot_wt);
	}
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container(this->bb_pot_set);
	if (bb_as_geom_pot_changed){
		//cout << "bb " << ca_as_geom_pot_wt << endl;
		bb_min_pot->get_default_energies_map().adjust_weight(potentials_name("as_geom_pot"), bb_as_geom_pot_wt);
	}

	LOOP::remodel_loops(pose_,
			this->resNum1,
			this->resNum2,
			//loop_mask,
			this->start_beta_steps,
			this->anneal_steps,
			this->final_beta_steps,
			this->start_beta,
			this->final_beta,
			100,
			ca_pot,//ca_pot_set,
			bb_min_pot);//bb_pot_set); // min_steps

	//LOOP::restore_backup_coords(pose_, stored_backup);

	double move_dist = 0;

	for (LOOP::atom_shared_ptr_vector3d_map::iterator i = stored_backup.begin(); i != stored_backup.end(); i++ ){
			const double dist = (i->second - i->first->get_coords()).mod();
			if (dist > move_dist){
				move_dist = dist;
			}
	}



	return_flags.move_completed = true;
	return_flags.move_dist = move_dist;
	return_flags.is_large_move = false;

	return true;
}


bb_composite_ca_mover::bb_composite_ca_mover(const int resNum1_,
		const int resNum2_,
		const unsigned long high_T_steps_,
		const unsigned long anneal_steps_,
		const unsigned long low_T_steps_,
		const double start_beta_,
		const double final_beta_) : resNum1(resNum1_),
		resNum2(resNum2_),
		start_beta_steps(high_T_steps_),
		anneal_steps(anneal_steps_),
		final_beta_steps(low_T_steps_),
		start_beta(start_beta_),
		final_beta(final_beta_),
		ca_pot_set("ca_default"),
		bb_pot_set("bb_min_default"){

	ca_as_geom_pot_changed =  false;
	bb_as_geom_pot_changed = false;

	ca_as_geom_pot_wt = 1.0;
	bb_as_geom_pot_wt = 1.0;

}

bb_composite_ca_mover_shared_ptr new_bb_composite_ca_mover(const int resNum1,
			const int resNum2,
			const unsigned long high_T_steps,
			const unsigned long anneal_steps,
			const unsigned long low_T_steps,
			const double start_beta,
			const double final_beta){

	if (resNum2 > resNum1){

		bb_composite_ca_mover_shared_ptr ptr(new bb_composite_ca_mover(resNum1, resNum2,
				high_T_steps,
				anneal_steps,
				low_T_steps,
				start_beta,
				final_beta));
		return ptr;

	}
	std::cerr << "ERROR: new_bb_composite_ca_mover: can not make mover\t" << resNum1 << "\t" << resNum2 << std::endl;
	return bb_composite_ca_mover_shared_ptr();
}

move_set_shared_ptr bb_composite_ca_mover_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const int min_seq_sep,
		const int max_seq_sep,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
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
				for (int i = resNum; i <= resNum + seq_sep; i++){
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

					/*
					std::cout << "made: " << resNum << "\t"
							<< resNum + seq_sep  << "\t"
							<< allowed_residues[resNum] << "\t"
							<< allowed_residues[resNum + seq_sep]
							<< std::endl;
					 */

					mover_shared_ptr ptr = new_bb_composite_ca_mover(resNum,
							resNum+seq_sep,
							high_T_steps,
							anneal_steps,
							low_T_steps,
							start_beta,
							final_beta);
					this_move_set->add_move(ptr, 1.0);
				}

			}


		}

	}


	return this_move_set;

}

move_set_shared_ptr bb_composite_ca_mover_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const int min_seq_sep,
		const int max_seq_sep,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta){

	pose_shared_ptr pose_ = meta_data->get_pose();
	bool_vector allowed_residues(pose_->get_residue_count(), true);

	return bb_composite_ca_mover_move_set_factory(meta_data,
			min_seq_sep,
			max_seq_sep,
			high_T_steps,
			anneal_steps,
			low_T_steps,
			start_beta,
			final_beta,
			allowed_residues);

}








}
}
}
}



