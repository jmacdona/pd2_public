/*
 * ca_radgyr.h
 *
 *  Created on: 8 Sep 2010
 *      Author: jmacdona
 */

#ifndef CA_RADGYR_H_
#define CA_RADGYR_H_

#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

potential_shared_ptr new_ca_radgyr();

class ca_radgyr : public potential_interface {

	friend potential_shared_ptr new_ca_radgyr();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

	//std::istream& load_data( std::istream& input );

	double mean_native_radgyr( const PRODART::POSE::const_pose_shared_ptr, const int chain) const;
	double stddev_native_radgyr( const PRODART::POSE::const_pose_shared_ptr, const int chain) const;

private:
	ca_radgyr();

	double mean_native_radgyr(const double resCount) const;
	double mean_random_walk_radgyr(const double resCount) const;

	double stddev_native_radgyr(const double resCount) const;
	double stddev_random_walk_radgyr(const double resCount) const;

	double k_native_radgyr(const double resCount) const;
	double k_random_walk_radgyr(const double resCount) const;

	//const PRODART::POSE::POTENTIALS::potentials_name name;

	double native_mean_a0, native_mean_a1;
	double rand_walk_mean_a0, rand_walk_mean_a1;

	double native_stddev_a0, native_stddev_a1;
	double rand_walk_stddev_a0, rand_walk_stddev_a1;


};

potential_shared_ptr new_ca_radgyr();

}
}
}
}

#endif /* CA_RADGYR_H_ */
