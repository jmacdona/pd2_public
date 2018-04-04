/*
 * residue_contact_restraint.cpp
 *
 *  Created on: Jul 4, 2012
 *      Author: jmacdona
 */


#include "residue_contact_restraint.h"

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

potential_shared_ptr new_residue_contact_restraint(){
	potential_shared_ptr ptr(new residue_contact_restraint());
	return ptr;
}

residue_contact_restraint::residue_contact_restraint() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("res_contact_rst"));
}






bool residue_contact_restraint::init(){

	return true;
}


double residue_contact_restraint::get_bond_energy(PRODART::POSE::META::upper_lower_bonded_pair_element& ele) const{
	if (ele.atom1_ptr->isActiveAndSet() && ele.atom2_ptr->isActiveAndSet() && (ele.atom1_ptr != ele.atom2_ptr)){
		const double dist = (ele.atom1_ptr->get_coords() -  ele.atom2_ptr->get_coords()).mod();
		/*
		cout << "TEST: " << ele.res_num1 << " " << ele.res_num2 << " " << dist << " "
				<< ele.equilibrium_dist_lower << " " << ele.equilibrium_dist_upper << " : ";
		 */
		if (dist < ele.equilibrium_dist_lower || dist > ele.equilibrium_dist_upper){
			const double equil_dist = dist <= ele.equilibrium_dist_lower ? ele.equilibrium_dist_lower : ele.equilibrium_dist_upper;
			const double energy = ele.half_bond_const * std::pow(dist - equil_dist , 2);
			ele.cached_score = energy;
			//cout << energy << "\n";
			return energy;
		}
		//cout << "\n";
	}
	return 0;
}

void residue_contact_restraint::bond_grad_mag(const double bond_length, double& energy, double& grad_mag, const double half_bond_const, const double equilibrium_dist) const{
	const double diff = bond_length - equilibrium_dist;
	grad_mag = (half_bond_const * 2.0 * (diff));
	energy = (half_bond_const * pow(diff, 2));
}

double residue_contact_restraint::bond_energy_grad(PRODART::POSE::META::upper_lower_bonded_pair_element& ele,
		PRODART::UTILS::vector3d_vector& grad,
		const double weight) const{

	double energy = 0, grad_mag = 0;
	PRODART::POSE::const_atom_shared_ptr atm1 = ele.atom1_ptr;
	PRODART::POSE::const_atom_shared_ptr atm2 = ele.atom2_ptr;

	if (atm1->isActiveAndSet() && atm2->isActiveAndSet() && (atm1 != atm2)){


		const double bond_len = (atm1->get_coords() - atm2->get_coords()).mod();
		if (bond_len < ele.equilibrium_dist_lower || bond_len > ele.equilibrium_dist_upper){
			const double equil_dist = bond_len <= ele.equilibrium_dist_lower ? ele.equilibrium_dist_lower : ele.equilibrium_dist_upper;
			bond_grad_mag(bond_len, energy, grad_mag, ele.half_bond_const, equil_dist);
			grad_mag *= weight;
			const PRODART::UTILS::vector3d unit_grad_vec = (atm1->get_coords() - atm2->get_coords()) / bond_len;
			const PRODART::UTILS::vector3d atm1_grad = grad_mag * unit_grad_vec;
			const PRODART::UTILS::vector3d atm2_grad = -grad_mag * unit_grad_vec;

			const int index1 = atm1->get_seq_num();
			const int index2 = atm2->get_seq_num();
			grad[index1] += atm1_grad;
			grad[index2] += atm2_grad;
		}

	}


	return energy;
}


double residue_contact_restraint::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::upper_lower_bonded_pair_element_vector& bond_list = pose_meta_->get_residue_contact_restraint_list();
	PRODART::POSE::META::upper_lower_bonded_pair_element_vector::iterator iter;

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

	//cout << "TEST: " << total_energy << "\n";

	return energies_map.add_energy_component(name_vector[0], total_energy);

}

double residue_contact_restraint::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);


	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::upper_lower_bonded_pair_element_vector& bond_list = pose_meta_->get_residue_contact_restraint_list();
	PRODART::POSE::META::upper_lower_bonded_pair_element_vector::iterator iter;

	double  total_energy = 0;

	for (iter = bond_list.begin(); iter != bond_list.end(); iter++){
		total_energy += this->bond_energy_grad(*iter, grad, weight);
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);

}


}
}
}


