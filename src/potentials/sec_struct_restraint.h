/*
 * sec_struct_restraint.h
 *
 *  Created on: 27 Jan 2011
 *      Author: jmacdona
 */

#ifndef SEC_STRUCT_RESTRAINT_H_
#define SEC_STRUCT_RESTRAINT_H_



#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_sec_struct_restraint();

class sec_struct_restraint : public potential_interface {


	friend potential_shared_ptr new_sec_struct_restraint();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

protected:
	sec_struct_restraint();



};



}
}
}




#endif /* SEC_STRUCT_RESTRAINT_H_ */
