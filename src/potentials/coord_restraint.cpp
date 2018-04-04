/*
 * coord_restraint.cpp
 *
 *  Created on: 3 Nov 2011
 *      Author: jmacdona
 */
#include "coord_restraint.h"

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

potential_shared_ptr new_coord_restraint(){
	potential_shared_ptr ptr(new coord_restraint());
	return ptr;
}

coord_restraint::coord_restraint() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("coord_restraint"));
}



double coord_restraint::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::coord_rst_element_vector& bond_list = pose_meta_->get_coord_restraint_list();
	PRODART::POSE::META::coord_rst_element_vector::iterator iter;

	double  total_energy = 0;

	for (iter = bond_list.begin(); iter != bond_list.end(); iter++){
		total_energy += this->get_bond_energy(*iter);
	}

	/* WRONG
	potentials_store& store = energies_map[base_name];
	store.energy = total_energy;
	energies_map.total_energy += store.weight * total_energy;
	return energies_map.total_energy;
	*/

	return energies_map.add_energy_component(name_vector[0], total_energy);

}

//TODO
double coord_restraint::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);


	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::coord_rst_element_vector& bond_list = pose_meta_->get_coord_restraint_list();
	PRODART::POSE::META::coord_rst_element_vector::iterator iter;

	double  total_energy = 0;

	for (iter = bond_list.begin(); iter != bond_list.end(); iter++){
		total_energy += this->bond_energy_grad(*iter, grad, weight);
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);

}


bool coord_restraint::init(){

	return true;
}


}
}
}


