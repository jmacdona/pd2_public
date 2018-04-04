/*
 * bb_old_hbond.h
 *
 *  Created on: 26 Oct 2010
 *      Author: jmacdona
 */

#ifndef BB_OLD_HBOND_H_
#define BB_OLD_HBOND_H_


#include "pose_meta/bb_pose_meta.h"
#include "../potential_interface.h"
#include <map>
#include <vector>

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{


typedef std::vector<double_vector> double_vector_vector;

potential_shared_ptr new_bb_old_hbond();

class bb_old_hbond : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_bb_old_hbond();

protected:

	enum hb_type {
		hb_other = 0,
				hb_helix = 1,
				hb_strand = 2,
				hb_not = 3
	};

	std::string hb_type_to_string(hb_type type) const{
		if (type == hb_other){
			return std::string("hb_other");
		}
		else if (type == hb_helix){
			return std::string("hb_helix");
		}
		else if (type == hb_strand){
			return std::string("hb_strand");
		}
		return std::string("hb_not");
	}

	static const int min_seq_sep = 3;

	static const double dist_bin_size;
	static const double psi_bin_size;
	static const double omega_bin_size;
	static const double tau_bin_size;

	static const double max_dist,
			min_dist;
	static const int max_dist_bin,
		min_dist_bin,
		max_psi_bin,
		max_omega_bin,
		max_tau_bin;

	static const double fade_off_dist;
	static const double angular_fade_off_dist;
	static const double zero_dist_energy;

	static const double grad_h;

	double_vector_vector H_O_dist_bin_score,
						  psi_bin_score,
						  omega_bin_score,
						  tau_bin_score;
	double_vector hbond_energy;

	int get_dist_bin(double) const;
	int get_omega_bin(double) const;
	int get_psi_bin(double) const;
	int get_tau_bin(double) const;

	double get_cdist_from_bin(int) const;
	double get_comega_from_bin(int) const;
	double get_cpsi_from_bin(int) const;
	double get_ctau_from_bin(int) const;

	double get_min_dist_from_bin(int) const;
	double get_min_omega_from_bin(int) const;
	double get_min_psi_from_bin(int) const;
	double get_min_tau_from_bin(int) const;

	bb_old_hbond();

	hb_type get_hbond_type(const PRODART::POSE::META::nb_pair_element& ele,
			const PRODART::POSE::META::bb_pose_meta_shared_ptr bb_meta_dat) const;

	double get_dist_energy(hb_type type, const double dist) const;
	double get_omega_energy(hb_type type, const double omega) const;
	double get_psi_energy(hb_type type, const double psi) const;
	double get_tau_energy(hb_type type, const double tau) const;

	double get_energy(const hb_type type, const double dist, const double omega, const double psi, const double tau) const;

	double get_energy(const hb_type type,
			const double dist,
			const UTILS::vector3d n,
			const UTILS::vector3d h,
			const UTILS::vector3d ca,
			const UTILS::vector3d c,
			const UTILS::vector3d o) const;

	double get_angle_fade_off(const double dist) const;

	double get_energy(const PRODART::POSE::META::nb_pair_element& ele,
			hb_type type,
			const PRODART::POSE::META::bb_pose_meta_shared_ptr bb_meta_dat) const;
	double get_energy_grad(const PRODART::POSE::META::nb_pair_element& ele,
			hb_type type,
			const PRODART::POSE::META::bb_pose_meta_shared_ptr bb_meta_dat,
			UTILS::vector3d_vector& grad,
			const double weight) const;


public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	std::istream& load_data( std::istream& input );

private:

	// force val into memory
	void dummy_funct(double val) const;



};






}
}
}
}


#endif /* BB_OLD_HBOND_H_ */
