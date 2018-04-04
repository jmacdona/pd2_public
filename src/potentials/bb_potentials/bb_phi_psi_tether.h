/*
 * bb_phi_psi_tether.h
 *
 *  Created on: 29 Oct 2010
 *      Author: jmacdona
 */

#ifndef BB_PHI_PSI_TETHER_H_
#define BB_PHI_PSI_TETHER_H_
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

potential_shared_ptr new_bb_phi_psi_tether();

class bb_phi_psi_tether : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_bb_phi_psi_tether();

protected:

	bb_phi_psi_tether();

	double get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const;
	double get_dE_ddih(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const;
	double get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele) const;
	double get_energy_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;
	double get_energy_ana_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;

	static const double h;// = 1e-8;

	double get_dihedral_distance(const double dih1, const double dih2) const;


	// force val into memory
	void dummy_funct(double val) const;

	mutable META::double_3_tuple_vector equil_vals;
	mutable bool equil_set;
	mutable META::simple_harmonic_dihedral_element_vector dih_vec;
	//mutable META::simple_harmonic_dihedral_element_vector psi_atoms;

	void set_equil(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_) const;

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	void reset_equil(){
		equil_set = false;
	}

	potential_shared_ptr get_new_instance() const{
		return new_bb_phi_psi_tether();
	}

private:



};

}
}
}
}


#endif /* BB_PHI_PSI_TETHER_H_ */
