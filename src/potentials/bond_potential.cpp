//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * bond_potential.cpp
 *
 *  Created on: 11 Jun 2010
 *      Author: jmacdona
 */

#include "bond_potential.h"

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

potential_shared_ptr new_bond_potential(){
	potential_shared_ptr ptr(new bond_potential());
	return ptr;
}

bond_potential::bond_potential() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bond_potential"));
}






/*
void bond_potential::bond_grad_mag(const double bond_length, double& energy, double& grad_mag, const double half_bond_const, const double equilibrium_dist) const{

		const double diff = bond_length - equilibrium_dist;
		grad_mag = (half_bond_const * 2.0 * (diff));
		energy = (half_bond_const * pow(diff, 2));

}
*/


/*
double bond_potential::bond_energy_grad(PRODART::POSE::META::bonded_pair_element& ele,
		UTILS::vector3d_vector& grad,
		const double weight) const{

	PRODART::POSE::const_atom_shared_ptr atm1 = ele.atom1_ptr;
	PRODART::POSE::const_atom_shared_ptr atm2 = ele.atom2_ptr;

	double energy = 0, grad_mag = 0;
	const double bond_len = (atm1->get_coords() - atm2->get_coords()).mod();
	bond_grad_mag(bond_len, energy, grad_mag, ele.half_bond_const, ele.equilibrium_dist);
	grad_mag *= weight;
	const vector3d unit_grad_vec = (atm1->get_coords() - atm2->get_coords()) / bond_len;
	const vector3d atm1_grad = grad_mag * unit_grad_vec;
	const vector3d atm2_grad = -grad_mag * unit_grad_vec;

	const int index1 = atm1->get_seq_num();
	const int index2 = atm2->get_seq_num();
	grad[index1] += atm1_grad;
	grad[index2] += atm2_grad;


	return energy;
}
*/


//! does it provide energies labeled with this name
/*
bool bond_potential::provides(const potentials_name& query_name) const{
	if (query_name == base_name){
		return true;
	}

	return false;
}
*/


bool bond_potential::init(){

	return true;
}


}
}
}


