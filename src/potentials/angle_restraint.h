/*
 * angle_restraint.h
 *
 *  Created on: 21 Feb 2011
 *      Author: jmacdona
 */

#ifndef ANGLE_RESTRAINT_H_
#define ANGLE_RESTRAINT_H_
#include "pose_meta/bb_pose_meta.h"
#include "potential_interface.h"
//#include <boost/tuple/tuple.hpp>
//#include <boost/tuple/tuple_comparison.hpp>
#include <map>

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

//typedef std::map<POSE::atom_type, double> atom_type_double_map;
//typedef boost::tuple<POSE::atom_type, POSE::atom_type> atom_type2_tuple;
//typedef std::map<atom_type2_tuple, double> atom_type2_tuple_double_map;

potential_shared_ptr new_angle_restraint();

class angle_restraint : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_angle_restraint();

protected:

	angle_restraint();

	double get_energy(const double ang, const double equilibrium_angle, const double half_angle_const) const;
	double get_energy(const PRODART::POSE::META::angle_element& ele) const;
	double get_energy_grad(const PRODART::POSE::META::angle_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;
	double get_dE_dang(const double ang, const double equilibrium_angle, const double half_angle_const) const;
	double get_energy_ana_grad(const PRODART::POSE::META::angle_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;

	static const double h;// = 1e-8;

	// force val into memory
	void dummy_funct(double val) const;

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

private:



};

}
}
}



#endif /* ANGLE_RESTRAINT_H_ */
