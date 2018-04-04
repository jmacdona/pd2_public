//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * potentials_container.cpp
 *
 *  Created on: 11 Jun 2010
 *      Author: jmacdona
 */

#include "potentials_container.h"

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


potentials_container_shared_ptr new_potentials_container(){
	potentials_container_shared_ptr ptr(new potentials_container());
	return ptr;
}

potentials_container::potentials_container(){

}


double potentials_container::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	assert(pose_meta_->get_pose()->coords_numerically_ok());

	//const const_pose_shared_ptr pose_ = _pose_meta->get_pose();
	potential_shared_ptr_vector::const_iterator iter;

	//double total_energy = 0;
	energies_map.total_energy = 0;

	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		//total_energy +=
		//assert(pose_->coords_numerically_ok());
		(*iter)->get_energy(pose_meta_, energies_map);
	}

	//assert(pose_meta_->get_pose()->coords_numerically_ok());

	return energies_map.total_energy;//total_energy;
}

double potentials_container::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	assert(pose_->coords_numerically_ok());

	/* do this when setting up pose_meta_
	vector3d_vector& grad = _pose_meta->get_gradient();
	grad.clear();
	grad.resize(pose_->get_all_atom_count(), vector3d(0,0,0));
	 */

	potential_shared_ptr_vector::const_iterator iter;

	vector3d_vector& grad = pose_meta_->get_gradient();
	vector3d_vector::iterator v3d_iter;
	const vector3d null_vec(0,0,0);
	//cout << "gradient size: " << grad.size() << endl;
	for (v3d_iter = grad.begin(); v3d_iter != grad.end(); v3d_iter++ ){
		*v3d_iter = null_vec;
	}
	//double total_energy = 0;
	energies_map.total_energy = 0;

	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		//total_energy +=
		//assert(pose_->coords_numerically_ok());
		(*iter)->get_energy_with_gradient(pose_meta_, energies_map);
	}

	assert(UTILS::is_vector3d_vector_valid(grad));

	//assert(pose_->coords_numerically_ok());

	return energies_map.total_energy;//total_energy;

}

void potentials_container::get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		bool_vector& vec) const{
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();

	potential_shared_ptr_vector::const_iterator iter;

	vec.clear();
	vec.resize(pose_->get_residue_count(), false);

	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		//total_energy +=
		(*iter)->get_residue_hot_spots(pose_meta_, vec);
	}

}

double potentials_container::get_energy_residue_loop(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map,
			const bool_vector& res_loop_mask) const{
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	assert(pose_meta_->get_pose()->coords_numerically_ok());

	//const const_pose_shared_ptr pose_ = _pose_meta->get_pose();
	potential_shared_ptr_vector::const_iterator iter;

	//double total_energy = 0;
	energies_map.total_energy = 0;

	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		//total_energy +=
		//assert(pose_->coords_numerically_ok());
		(*iter)->get_energy_residue_loop(pose_meta_, energies_map, res_loop_mask);
	}

	//assert(pose_meta_->get_pose()->coords_numerically_ok());

	return energies_map.total_energy;//total_energy;
}


bool potentials_container::init(){

	potential_shared_ptr_vector::const_iterator iter;

	bool rtn = true;

	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		rtn = rtn && (*iter)->init();

	}

	return rtn;

}


//! does it provide energies labeled with this name
bool potentials_container::provides(const potentials_name& query_name) const{

	potential_shared_ptr_vector::const_iterator iter;

	bool rtn = false;;

	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		if ( (*iter)->provides(query_name) == true){
			return true;
		}
	}

	return rtn;

}


void potentials_container::add_potential(potential_shared_ptr ptr){
	if (!ptr) return;
	potential_shared_ptr_vector::const_iterator iter;
	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		// check not already added
		if (*iter == ptr){
			return;
		}
	}
	if (ptr) potentials.push_back(ptr);
}

void potentials_container::add_potential(potential_shared_ptr ptr, const double weight){
	if (!ptr) return;
	potential_shared_ptr_vector::const_iterator iter;
	bool already_in = false;
	for (iter = potentials.begin(); iter != potentials.end(); iter++){
		// check not already added
		if (*iter == ptr){
			already_in = true;
		}
	}
	if (!already_in){
		if (ptr) potentials.push_back(ptr);
	}

	potentials_name_vector names = ptr->provides();
	potentials_name_vector::const_iterator nit;
	for (nit = names.begin(); nit != names.end(); nit++){
		this->default_energies_map[*nit].energy = 0;
		this->default_energies_map[*nit].weight = weight;
	}


}


potential_shared_ptr potentials_container::get_potential(const potentials_name& name){
	potential_shared_ptr_vector::iterator pot_iter;
	for (pot_iter = potentials.begin(); pot_iter != potentials.end(); pot_iter++){
		if ((*pot_iter)->provides(name)){
			return (*pot_iter);
		}
	}
	return potential_shared_ptr();
}


}
}
}

