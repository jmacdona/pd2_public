/*
 * bb_frag3_mq.h
 *
 *  Created on: Oct 31, 2010
 *      Author: jmacdon
 */

#ifndef BB_FRAG3_MQ_H_
#define BB_FRAG3_MQ_H_

#include "pose_meta/bb_pose_meta.h"
#include "../potential_interface.h"
#include <map>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

typedef std::map< std::string, double, std::less< std::string > > string_double_map;

potential_shared_ptr new_bb_frag3_mq();

class bb_frag3_mq : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_bb_frag3_mq();

protected:

	bb_frag3_mq();

	int get_dih_bin(const double dih, const int num_bins) const;

	string_double_map frag3_counts;

	int num_phi_psi_bins;
	int num_omega_bins;
	int count_cutoff;



public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	void get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			bool_vector& vec) const;

	double get_energy_residue_loop(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
				potentials_energies_map& energies_map,
				const bool_vector& res_loop_mask) const;

	bool init();

	std::istream& load_data( std::istream& input );

private:



};

}
}
}
}



#endif /* BB_FRAG3_MQ_H_ */
