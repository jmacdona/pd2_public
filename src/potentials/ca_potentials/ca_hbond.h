/*
 * ca_hbond.h
 *
 *  Created on: 1 Sep 2010
 *      Author: jmacdona
 */

#ifndef CA_HBOND_H_
#define CA_HBOND_H_
#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

typedef std::vector< double > double_vector;
typedef std::vector< double_vector > double_vector_vector;

potential_shared_ptr new_ca_hbond();

class ca_hbond : public potential_interface {


	friend potential_shared_ptr new_ca_hbond();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

	std::istream& load_data( std::istream& input );

private:
	//const PRODART::POSE::POTENTIALS::potentials_name name;
	//potentials_name_vector name_vector;
	ca_hbond();

	int getHbond_type(PRODART::POSE::four_state_sec_struct phipsiombin_i, PRODART::POSE::four_state_sec_struct phipsiombin_j , int seq_sep) const;
	int get_N_O_dist_bin(double) const;
	int get_chi_i_bin(double) const;
	int get_chi_j_bin(double) const;
	int get_tau_bin(double) const;

	double_vector_vector N_O_dist_bin_score,
		chi_i_bin_score,
		chi_j_bin_score,
		tau_bin_score;

	double_vector hbond_energy;

};












}
}
}
}

#endif /* CA_HBOND_H_ */

