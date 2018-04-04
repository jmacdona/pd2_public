//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * ca_bump.h
 *
 *  Created on: Jun 6, 2010
 *      Author: jmacdon
 */

#ifndef CA_BUMP_H_
#define CA_BUMP_H_

#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

potential_shared_ptr new_ca_bump();

class ca_bump : public potential_interface {


	friend potential_shared_ptr new_ca_bump();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

private:
	ca_bump();

	double get_energy(const PRODART::POSE::META::nb_pair_element& ele) const;
	double get_energy_grad(const PRODART::POSE::META::nb_pair_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;


	double CA_radius;
	double CA_radius_sq;

	//const PRODART::POSE::POTENTIALS::potentials_name name;

};


}
}
}
}

#endif /* CA_BUMP_H_ */
