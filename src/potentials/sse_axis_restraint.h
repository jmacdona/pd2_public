/*
 * sse_axis_restraint.h
 *
 *  Created on: Mar 13, 2012
 *      Author: jmacdona
 */

#ifndef SSE_AXIS_RESTRAINT_H_
#define SSE_AXIS_RESTRAINT_H_

#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"

#include "pose_utils/sse_axes.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_sse_axis_restraint();

class sse_axis_restraint : public potential_interface {


	friend potential_shared_ptr new_sse_axis_restraint();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;


	double get_energy(const META::sse_axis_element & ele, UTILS::vector3d_vector ca_coords) const;


protected:
	sse_axis_restraint();

	// force val into memory
	void dummy_funct(double val) const;
	static const double h;




};

}
}
}
#endif /* SSE_AXIS_RESTRAINT_H_ */
