/*
 * ca_bond.h
 *
 *  Created on: 26 Aug 2010
 *      Author: jmacdona
 */

#ifndef CA_BOND_H_
#define CA_BOND_H_


#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

potential_shared_ptr new_ca_bond();

class ca_bond : public potential_interface {


	friend potential_shared_ptr new_ca_bond();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

	std::istream& load_data( std::istream& input );

private:
	ca_bond();
	//const PRODART::POSE::POTENTIALS::potentials_name name;

	double bond_energy(PRODART::POSE::four_state_sec_struct, const double bond_length) const;
	void bond_grad_mag(PRODART::POSE::four_state_sec_struct, const double bond_length, double& energy, double& grad_mag) const;

	double bond_energy_grad(PRODART::POSE::four_state_sec_struct,
			PRODART::POSE::const_atom_shared_ptr atm1,
			PRODART::POSE::const_atom_shared_ptr atm2,
			UTILS::vector3d_vector& grad,
			const double weight) const;

	double frag_bond_energy(const PRODART::POSE::META::frag4_element& ele) const;
	double frag_bond_energy_grad(const PRODART::POSE::META::frag4_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;

    double trans_CA_CA_l0;
    double cis_CA_CA_l0;

    double trans_CA_CA_k;
    double cis_CA_CA_k;


};


}
}
}
}





#endif /* CA_BOND_H_ */
