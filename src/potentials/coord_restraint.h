/*
 * coord_restraint.h
 *
 *  Created on: 3 Nov 2011
 *      Author: jmacdona
 */

#ifndef COORD_RESTRAINT_H_
#define COORD_RESTRAINT_H_
#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_coord_restraint();

class coord_restraint : public potential_interface {


	friend potential_shared_ptr new_coord_restraint();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

protected:
	coord_restraint();


	//PRODART::POSE::POTENTIALS::potentials_name base_name;

	double get_bond_energy(PRODART::POSE::META::coord_rst_element& ele) const{
		if (ele.atom_ptr->isActiveAndSet()){
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

	double bond_energy_grad(PRODART::POSE::META::coord_rst_element& ele,
			PRODART::UTILS::vector3d_vector& grad,
			const double weight) const{
		const double tol = 0.0;
		double energy = 0, grad_mag = 0;
		PRODART::POSE::const_atom_shared_ptr atm1 = ele.atom_ptr;


		if (atm1->isActiveAndSet()){


			const double bond_len = ele.dist;//(atm1->get_coords() - atm2->get_coords()).mod();
			if (fabs(ele.dist) > tol){
				bond_grad_mag(bond_len, energy, grad_mag, ele.half_bond_const, ele.equilibrium_dist);
				grad_mag *= weight;
				const PRODART::UTILS::vector3d unit_grad_vec = (atm1->get_coords() - ele.coord) / bond_len;
				const PRODART::UTILS::vector3d atm1_grad = grad_mag * unit_grad_vec;


				const int index1 = atm1->get_seq_num();

				grad[index1] += atm1_grad;
			}


		}


		return energy;
	}



};







}
}
}

#endif /* COORD_RESTRAINT_H_ */
