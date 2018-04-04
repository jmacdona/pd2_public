/*
 * bb_dihedral.h
 *
 *  Created on: 27 Oct 2010
 *      Author: jmacdona
 */

#ifndef BB_DIHEDRAL_H_
#define BB_DIHEDRAL_H_
#include "pose_meta/bb_pose_meta.h"
#include "../potential_interface.h"
//#include <boost/tuple/tuple.hpp>
//#include <boost/tuple/tuple_comparison.hpp>
#include <map>
#include "utils/dihedral_derivative.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

//typedef std::map<POSE::atom_type, double> atom_type_double_map;
//typedef boost::tuple<POSE::atom_type, POSE::atom_type> atom_type2_tuple;
//typedef std::map<atom_type2_tuple, double> atom_type2_tuple_double_map;

potential_shared_ptr new_bb_dihedral();

class bb_dihedral : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_bb_dihedral();

protected:

	bb_dihedral();

	double get_energy(const PRODART::POSE::META::dihedral_element& ele, const double dih) const;
	double get_dE_ddih(const PRODART::POSE::META::dihedral_element& ele, const double dih) const;
	double get_energy(const PRODART::POSE::META::dihedral_element& ele) const;
	double get_dE_ddih(const PRODART::POSE::META::dihedral_element& ele) const;
	double get_energy_grad(const PRODART::POSE::META::dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;
	double get_energy_ana_grad(const PRODART::POSE::META::dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;

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
}


#endif /* BB_DIHEDRAL_H_ */
