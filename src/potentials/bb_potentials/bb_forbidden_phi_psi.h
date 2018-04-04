/*
 * bb_forbidden_phi_psi.h
 *
 *  Created on: 1 Nov 2010
 *      Author: jmacdona
 */

#ifndef BB_FORBIDDEN_PHI_PSI_H_
#define BB_FORBIDDEN_PHI_PSI_H_

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


potential_shared_ptr new_bb_forbidden_phi_psi();

class bb_forbidden_phi_psi : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_bb_forbidden_phi_psi();

protected:

	bb_forbidden_phi_psi();

	void getPhiPsiBin(const int phi_psi_sector, int& phi, int& psi) const;
	double get_lower_bound(int bin) const;
	double get_grid_centre(int bin) const;

	double_vector overall_freq_map;
	double_vector energies;

	double totalResidueCount;

	static const unsigned int grid_divs;

	double phi_psi_increment;

	double forbidden_cutoff;
	double hot_spot_cutoff;

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




#endif /* BB_FORBIDDEN_PHI_PSI_H_ */
