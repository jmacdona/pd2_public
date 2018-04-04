/*
 * bb_bond.cpp
 *
 *  Created on: 25 Oct 2010
 *      Author: jmacdona
 */
#include "bb_bond.h"

using namespace boost;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;
using namespace std;

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

potential_shared_ptr new_bb_bond(){
	potential_shared_ptr ptr(new bb_bond());
	return ptr;
}

bb_bond::bb_bond(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_bond"));
}

double bb_bond::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::bonded_ele_vector& bond_list = bb_meta_dat->get_bb_bond_list();//pose_meta_->get_bond_list();
	PRODART::POSE::META::bonded_ele_vector::iterator iter;

	double  total_energy = 0;

	for (iter = bond_list.begin(); iter != bond_list.end(); iter++){
		total_energy += this->get_bond_energy(*iter);
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_bond::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);


	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::bonded_ele_vector& bond_list = bb_meta_dat->get_bb_bond_list();//pose_meta_->get_bond_list();
	PRODART::POSE::META::bonded_ele_vector::iterator iter;

	double  total_energy = 0;

	for (iter = bond_list.begin(); iter != bond_list.end(); iter++){
		total_energy += this->bond_energy_grad(*iter, grad, weight);
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}


bool bb_bond::init(){

	return true;
}




}
}
}
}



