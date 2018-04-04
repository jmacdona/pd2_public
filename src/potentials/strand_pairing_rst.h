/*
 * strand_pairing.h
 *
 *  Created on: 28 Jan 2011
 *      Author: jmacdona
 */

#ifndef STRAND_PAIRING_H_
#define STRAND_PAIRING_H_





#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_strand_pairing_rst();

// TODO try adding CB-CB pair dist restraint 3.5 to 6.2 angstr
class strand_pairing_rst : public potential_interface {


	friend potential_shared_ptr new_strand_pairing_rst();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

protected:
	strand_pairing_rst();

	double get_ca_pair_energy(POSE::pose_shared_ptr pose_,
			const int res1, const int res2,
			const bool is_parallel, const double weight) const;
	double get_cb_pair_energy(POSE::pose_shared_ptr pose_,
			const int res1, const int res2,
			const bool is_parallel, const double weight) const;
	double get_bb_pair_energy(POSE::pose_shared_ptr pose_,
			const int res1, const int res2,
			const bool is_parallel, const double weight) const;

	double get_ca_pair_energy_grad(POSE::pose_shared_ptr pose_,
			const int res1, const int res2,
			const bool is_parallel, const double weight,
			PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const;
	double get_cb_pair_energy_grad(POSE::pose_shared_ptr pose_,
			const int res1, const int res2,
			const bool is_parallel, const double weight,
			PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const;
	double get_bb_pair_energy_grad(POSE::pose_shared_ptr pose_,
			const int res1, const int res2,
			const bool is_parallel, const double weight,
			PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const;


	double get_dist_energy(const double lower_cutoff, const double upper_cutoff, const double dist) const;
	double get_dE_ddist(const double lower_cutoff, const double upper_cutoff, const double dist) const;
	double get_pair_energy_grad(const double lower_cutoff,
			const double upper_cutoff, POSE::atom_shared_ptr, POSE::atom_shared_ptr,
			const double weight,
			PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const;
	double get_pair_energy(const double lower_cutoff,
			const double upper_cutoff, POSE::atom_shared_ptr, POSE::atom_shared_ptr,
			const double weight) const;



};



}
}
}
























#endif /* STRAND_PAIRING_H_ */
