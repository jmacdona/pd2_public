//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * ca_bump.cpp
 *
 *  Created on: Jun 6, 2010
 *      Author: jmacdon
 */

#include "ca_bump.h"

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
namespace CA{


potential_shared_ptr new_ca_bump(){
	potential_shared_ptr ptr(new ca_bump());
	return ptr;
}

ca_bump::ca_bump() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_bump"));
	CA_radius = 4.0;
	CA_radius_sq = CA_radius * CA_radius;
}


double ca_bump::get_energy(const PRODART::POSE::META::nb_pair_element& ele) const{
	const double dist_sq = ele.dist_sq;
	//TODO sort out seq sep in meta
	const int seq_sep = ele.seq_sep;	//99999;
	if (dist_sq < CA_radius_sq
			&& seq_sep > 1){
		return (CA_radius_sq - (dist_sq)) / CA_radius;
		//total_score += pow((CA_radius_sq - (dist_sq)),2) / CA_radius;
	}
	return 0;
}

double ca_bump::get_energy_grad(const PRODART::POSE::META::nb_pair_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	const double dist_sq = ele.dist_sq;
	//TODO sort out seq sep in meta
	const int seq_sep = ele.seq_sep;	//99999;
	if (dist_sq < CA_radius_sq
			&& seq_sep > 1){
		const double energy = (CA_radius_sq - (dist_sq)) / CA_radius;
		const double grad_mag = weight * ((2.0 * ele.dist) / CA_radius);

		const vector3d unit_grad_vec = (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()) / ele.dist;
		const vector3d atm1_grad = -grad_mag * unit_grad_vec;
		const vector3d atm2_grad = grad_mag * unit_grad_vec;

		const int index1 = ele.atom1_ptr->get_seq_num();
		const int index2 = ele.atom2_ptr->get_seq_num();

		grad[index1] += atm1_grad;
		grad[index2] += atm2_grad;

		return energy;
	}
	else{
		return 0;
	}
}

double ca_bump::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	double total_score = 0;

	const nb_ele_vector& pair_list = ca_meta_dat->get_ca_pair_list();
	nb_ele_vector::const_iterator iter;

	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){

		//TODO sort out seq sep in meta
		//const int seq_sep = iter->seq_sep;	//99999;
		total_score += this->get_energy(*iter);
		/*
		if (dist_sq < CA_radius_sq
				&& seq_sep > 1){
			//total_score += (CA_radius_sq - (dist_sq)) / CA_radius;
			//total_score += pow((CA_radius_sq - (dist_sq)),2) / CA_radius;
		}
		*/

	}

	/* WRONG
	potentials_store& store = energies_map[name];
	store.energy = total_score;
	energies_map.total_energy += store.weight * total_score;
	return energies_map.total_energy;
	*/

	return energies_map.add_energy_component(this->name_vector[0], total_score);
}

//TODO
double ca_bump::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	vector3d_vector& grad = _pose_meta->get_gradient();
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	double total_score = 0;
	const double weight = energies_map.get_weight(this->name_vector[0]);
	const nb_ele_vector& pair_list = ca_meta_dat->get_ca_pair_list();
	nb_ele_vector::const_iterator iter;

	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){

		total_score += this->get_energy_grad(*iter, grad, weight);


	}



	return energies_map.add_energy_component(this->name_vector[0], total_score);
}


bool ca_bump::init(){

	return true;
}

//! does it provide energies labeled with this name
/*
bool ca_bump::provides(const potentials_name& query_name) const{
	if (query_name == name){
		return true;
	}

	return false;
}
*/










}
}
}
}
