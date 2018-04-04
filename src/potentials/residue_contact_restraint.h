/*
 * residue_contact_restraint.h
 *
 *  Created on: Jul 4, 2012
 *      Author: jmacdona
 */

#ifndef RESIDUE_CONTACT_RESTRAINT_H_
#define RESIDUE_CONTACT_RESTRAINT_H_


#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_residue_contact_restraint();

class residue_contact_restraint : public potential_interface {



	friend potential_shared_ptr new_residue_contact_restraint();

public:

	virtual ~residue_contact_restraint(){
	}

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

protected:
	residue_contact_restraint();


	//PRODART::POSE::POTENTIALS::potentials_name base_name;

	double get_bond_energy(PRODART::POSE::META::upper_lower_bonded_pair_element& ele) const;

	void bond_grad_mag(const double bond_length, double& energy, double& grad_mag, const double half_bond_const, const double equilibrium_dist) const;

	double bond_energy_grad(PRODART::POSE::META::upper_lower_bonded_pair_element& ele,
			PRODART::UTILS::vector3d_vector& grad,
			const double weight) const;

};




}
}
}



#endif /* RESIDUE_CONTACT_RESTRAINT_H_ */
