/*
 * sec_struct_dih_ang_rst.h
 *
 *  Created on: 8 Nov 2011
 *      Author: jmacdona
 */

#ifndef SEC_STRUCT_DIH_ANG_RST_H_
#define SEC_STRUCT_DIH_ANG_RST_H_


#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"
#include "utils/dihedral_derivative.h"
#include "utils/angle_derivative.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_sec_struct_dih_ang_rst();

class sec_struct_dih_ang_rst : public potential_interface {


	friend potential_shared_ptr new_sec_struct_dih_ang_rst();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

protected:
	sec_struct_dih_ang_rst();

	double get_dih_energy(const double equil_dih_angle, const double dih) const;
	double get_dE_ddih(const double equil_dih_angle, const double dih) const;
	double get_dih_energy(POSE::atom_shared_ptr atom1_ptr,
			POSE::atom_shared_ptr atom2_ptr,
			POSE::atom_shared_ptr atom3_ptr,
			POSE::atom_shared_ptr atom4_ptr,
			const double dih,
			const double equil_ang) const;
	double get_dih_energy_ana_grad(POSE::atom_shared_ptr atom1_ptr,
			POSE::atom_shared_ptr atom2_ptr,
			POSE::atom_shared_ptr atom3_ptr,
			POSE::atom_shared_ptr atom4_ptr,
			double dih,
			const double equil_ang,
			UTILS::vector3d_vector& grad, const double weight ) const;

	double get_min_max_dih_energy(POSE::atom_shared_ptr atom1_ptr,
			POSE::atom_shared_ptr atom2_ptr,
			POSE::atom_shared_ptr atom3_ptr,
			POSE::atom_shared_ptr atom4_ptr,
			const double dih,
			const double min_dih,
			const double max_dih) const;
	double get_min_max_dih_energy_ana_grad(POSE::atom_shared_ptr atom1_ptr,
			POSE::atom_shared_ptr atom2_ptr,
			POSE::atom_shared_ptr atom3_ptr,
			POSE::atom_shared_ptr atom4_ptr,
			double dih,
			const double min_dih,
			const double max_dih,
			UTILS::vector3d_vector& grad, const double weight ) const;


	static const double beta_min_phi, beta_max_phi;
	static const double beta_min_psi, beta_max_psi;

	static const double alpha_min_phi, alpha_max_phi;
	static const double alpha_min_psi, alpha_max_psi;

	static const double beta_min_ca_dih, beta_max_ca_dih;
	static const double alpha_min_ca_dih, alpha_max_ca_dih;

};



}
}
}





#endif /* SEC_STRUCT_DIH_ANG_RST_H_ */
