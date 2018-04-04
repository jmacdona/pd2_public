//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * bond_potential.h
 *
 *  Created on: 11 Jun 2010
 *      Author: jmacdona
 */

#ifndef BOND_POTENTIAL_H_
#define BOND_POTENTIAL_H_


#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_bond_potential();

class bond_potential : public potential_interface {



	friend potential_shared_ptr new_bond_potential();

public:

	virtual ~bond_potential(){
	}

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

protected:
	bond_potential();


	//PRODART::POSE::POTENTIALS::potentials_name base_name;

	double get_bond_energy(PRODART::POSE::META::bonded_pair_element& ele) const{
		if (ele.atom1_ptr->isActiveAndSet() && ele.atom2_ptr->isActiveAndSet() && (ele.atom1_ptr != ele.atom2_ptr)){
			const double dist = ele.dist;//(ele.atom1_ptr->get_coords() -  ele.atom2_ptr->get_coords()).mod();
			const double energy = ele.half_bond_const * std::pow(dist - ele.equilibrium_dist , 2);
			ele.cached_score = energy;
			return energy;
		}
		return 0;
	}

	void bond_grad_mag(const double bond_length, double& energy, double& grad_mag, const double half_bond_const, const double equilibrium_dist) const{
		const double diff = bond_length - equilibrium_dist;
		grad_mag = (half_bond_const * 2.0 * (diff));
		energy = (half_bond_const * pow(diff, 2));
	}

	double bond_energy_grad(PRODART::POSE::META::bonded_pair_element& ele,
			PRODART::UTILS::vector3d_vector& grad,
			const double weight) const{

		double energy = 0, grad_mag = 0;
		PRODART::POSE::const_atom_shared_ptr atm1 = ele.atom1_ptr;
		PRODART::POSE::const_atom_shared_ptr atm2 = ele.atom2_ptr;

		if (atm1->isActiveAndSet() && atm2->isActiveAndSet() && (atm1 != atm2)){


			const double bond_len = ele.dist;//(atm1->get_coords() - atm2->get_coords()).mod();
			bond_grad_mag(bond_len, energy, grad_mag, ele.half_bond_const, ele.equilibrium_dist);
			grad_mag *= weight;
			const PRODART::UTILS::vector3d unit_grad_vec = (atm1->get_coords() - atm2->get_coords()) / bond_len;
			const PRODART::UTILS::vector3d atm1_grad = grad_mag * unit_grad_vec;
			const PRODART::UTILS::vector3d atm2_grad = -grad_mag * unit_grad_vec;

			const int index1 = atm1->get_seq_num();
			const int index2 = atm2->get_seq_num();
			grad[index1] += atm1_grad;
			grad[index2] += atm2_grad;

		}


		return energy;
	}

};


inline double bond_potential::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::bonded_ele_vector& bond_list = pose_meta_->get_dist_harmonic_restraint_list();
	PRODART::POSE::META::bonded_ele_vector::iterator iter;

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
inline double bond_potential::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);


	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::bonded_ele_vector& bond_list = pose_meta_->get_dist_harmonic_restraint_list();
	PRODART::POSE::META::bonded_ele_vector::iterator iter;

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



#endif /* BOND_POTENTIAL_H_ */
