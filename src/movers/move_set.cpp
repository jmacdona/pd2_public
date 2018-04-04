//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * move_set.cpp
 *
 *  Created on: Mar 1, 2010
 *      Author: jmacdon
 */
#include "move_set.h"




namespace PRODART {
namespace POSE {
namespace MOVERS {

move_set_shared_ptr new_move_set(){
	boost::shared_ptr<move_set> nms(new move_set);
	//npose->_this = nms;
	return nms;
}

int move_set::max_tries = 100;

bool move_set::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{

	mover_flags& return_flags = meta_data->get_mover_flags();
	const PRODART::POSE::pose_shared_ptr protein = meta_data->get_pose();
	const unsigned int num_moves = this->move_vec.size();

	if (num_moves == 0){
		std::cout << "ERROR: move_set: no moves in set\n";
		return_flags.move_completed = false;
		return_flags.is_large_move = false;
		return false;
	}


	bool accepted = false;
	int tries = 0;
	while (accepted == false && tries < max_tries){
		unsigned int index = this->rand_gen->randInt(num_moves - 1);
		const double prob_acc = move_vec[index].prob_accept;
		if (prob_acc >= 1.0){
			return move_vec[index].mover->make_move(meta_data);
			accepted = true;
		}
		else if (rand_gen->rand() < prob_acc) {
			return move_vec[index].mover->make_move(meta_data);
			accepted = true;
		}
		tries++;
	}

	std::cout << "WARNING: move_set: exceeded tries limit\n";

	return_flags.move_completed = true;
	return_flags.is_large_move = false;

	return accepted;
}

void move_set::add_move(mover_shared_ptr ptr, const double weight){

	if (weight <= 0 ){
		throw move_set_bad_move_add();
	}

	move_set_entry new_ms_entry;
	new_ms_entry.mover = ptr;
	new_ms_entry.weight = weight;
	//new_ms_entry.prob_accept = 1.0;
	this->move_vec.push_back(new_ms_entry);
	this->calc_prob_accepts();
}

double move_set::find_max_weight() const{

	move_set_entry_vector::const_iterator iter;
	double max = -1;
	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		if (iter->weight > max){
			max =iter->weight;
		}
	}

	return max;
}

void move_set::calc_prob_accepts(){
	move_set_entry_vector::iterator iter;
	const double max = this->find_max_weight();
	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->prob_accept = iter->weight / max;
	}
}

void move_set::propagate_rand_num_gen(const MTRand::MTRand_shared_ptr ptr){
	move_set_entry_vector::iterator iter;

	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->mover->propagate_rand_num_gen(ptr);
	}
	rand_gen = ptr;
}

void move_set::propagate_start_beta(const double val){
	move_set_entry_vector::iterator iter;

	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->mover->propagate_start_beta(val);
	}
}
void move_set::propagate_final_beta(const double val){
	move_set_entry_vector::iterator iter;

	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->mover->propagate_final_beta(val);
	}
}
void move_set::propagate_ca_potential_preset(const std::string str){
	move_set_entry_vector::iterator iter;

	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->mover->propagate_ca_potential_preset(str);
	}
}
void move_set::propagate_bb_potential_preset(const std::string str){
	move_set_entry_vector::iterator iter;

	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->mover->propagate_bb_potential_preset(str);
	}
}
void move_set::propagate_ca_potential_weight(const POTENTIALS::potentials_name& name, const double weight){
	move_set_entry_vector::iterator iter;

	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->mover->propagate_ca_potential_weight(name, weight);
	}
}
void move_set::propagate_bb_potential_weight(const POTENTIALS::potentials_name& name, const double weight){
	move_set_entry_vector::iterator iter;

	for (iter = move_vec.begin(); iter != move_vec.end(); iter++){
		iter->mover->propagate_bb_potential_weight(name, weight);
	}
}

move_set::move_set(){
	this->type = mover_interface::mt_move_set;
	move_vec.reserve(1000);

	this->auto_rand_num_assign();
}

unsigned int move_set::get_move_count() const{
	return move_vec.size();
}


}
}
}

