/*
 * bb_strict_omega.h
 *
 *  Created on: 4 Nov 2011
 *      Author: jmacdona
 */

#ifndef BB_STRICT_OMEGA_H_
#define BB_STRICT_OMEGA_H_
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

potential_shared_ptr new_bb_strict_omega();


class bb_strict_omega : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_bb_strict_omega();

protected:

	bb_strict_omega();

	double get_energy(const double equil_dih_angle, const double dih) const;
	double get_dE_ddih(const double equil_dih_angle, const double dih) const;
	double get_energy(POSE::atom_shared_ptr atom1_ptr,
			POSE::atom_shared_ptr atom2_ptr,
			POSE::atom_shared_ptr atom3_ptr,
			POSE::atom_shared_ptr atom4_ptr,
			const double dih,
			const double equil_ang) const;
	double get_energy_ana_grad(POSE::atom_shared_ptr atom1_ptr,
			POSE::atom_shared_ptr atom2_ptr,
			POSE::atom_shared_ptr atom3_ptr,
			POSE::atom_shared_ptr atom4_ptr,
			double dih,
			const double equil_ang,
			UTILS::vector3d_vector& grad, const double weight ) const;

	//static const double h = 1e-8;

	double get_dihedral_distance(const double dih1, const double dih2) const;


	// force val into memory
	//void dummy_funct(double val) const;



	void set_equil(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_) const;

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

#endif /* BB_STRICT_OMEGA_H_ */
