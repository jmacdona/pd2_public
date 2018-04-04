/*
 * bb_pp_bias_correction.h
 *
 *  Created on: 1 Nov 2010
 *      Author: jmacdona
 */

#ifndef BB_PP_BIAS_CORRECTION_H_
#define BB_PP_BIAS_CORRECTION_H_

#include "pose_meta/bb_pose_meta.h"
#include "../potential_interface.h"
#include <map>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include "utils/hist2d_interp.h"
#include "utils/dihedral_derivative.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

class bb_pp_bias_correction;

typedef boost::shared_ptr<bb_pp_bias_correction> bb_pp_bias_correction_shared_ptr;

bb_pp_bias_correction_shared_ptr new_bb_pp_bias_correction();

class bb_pp_bias_correction : public POTENTIALS::potential_interface {

	friend bb_pp_bias_correction_shared_ptr new_bb_pp_bias_correction();

protected:

	bb_pp_bias_correction();

	void getPhiPsiBin(const int phi_psi_sector, int& phi, int& psi) const;
	double get_lower_bound(int bin) const;
	double get_grid_centre(int bin) const;


	double_vector overall_freq_map;
	double_vector energies;
	double totalResidueCount;
	static const unsigned int grid_divs;
	double phi_psi_increment;
	double hot_spot_cutoff;

	hist2d_interp histo;

	/*
	double_vector overall_freq_map;
	double_vector energies;
	double totalResidueCount;
	static const int grid_divs;
	double phi_psi_increment;
	double forbidden_cutoff;
	*/

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	void get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			bool_vector& vec) const;

	bool init();

	std::istream& load_data( std::istream& input );


	bool addPseudoCounts( const double sumCount = 100 );

	bool calculateScores( void );

	void addToDB(const double phi, const double psi);

    std::ostream& output_phi_psi_info(std::ostream& output);

private:





};

}
}
}
}





#endif /* BB_PP_BIAS_CORRECTION_H_ */
