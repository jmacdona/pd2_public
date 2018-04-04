/*
 * ca_density.h
 *
 *  Created on: 21 Oct 2013
 *      Author: jmacdona
 */

#ifndef CA_DENSITY_H_
#define CA_DENSITY_H_

#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"

#include <boost/tuple/tuple.hpp>


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{


potential_shared_ptr new_ca_density();

class ca_density : public potential_interface {


	friend potential_shared_ptr new_ca_density();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool training_pose(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_);
	bool background_training_pose(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_);
	bool calc_potentials_using_training();
	bool output_training_results(std::ostream& output);


	std::istream& load_data( std::istream& input );
	bool init();



private:

	double cutoff;
	int min_seq_sep;
	int max_density;

	std::map<int, double> alpha_counts, beta_counts, other_counts, all_counts;
	double total_alpha, total_beta, total_other, total_all;
	std::map<int, double> bg_alpha_counts, bg_beta_counts, bg_other_counts,  bg_all_counts;
	double bg_total_alpha, bg_total_beta, bg_total_other, bg_total_all;

	std::map<int, double> alpha_pot, beta_pot, other_pot, all_pot;

	ca_density();

	int_vector get_ca_counts(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_, const double cutoff) const;

};



}
}
}
}
#endif /* CA_DENSITY_H_ */
