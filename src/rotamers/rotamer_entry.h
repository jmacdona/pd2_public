/*
 * rotamer_entry.h
 *
 *  Created on: 28 Aug 2011
 *      Author: jmacdona
 */

#ifndef ROTAMER_ENTRY_H_
#define ROTAMER_ENTRY_H_

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "pose/residue_type.h"
#include "utils/vector3d.h"


namespace PRODART {
namespace ROTAMERS {

class rotamer_entry;

typedef boost::shared_ptr<rotamer_entry> rotamer_entry_shared_ptr;
typedef boost::shared_ptr<const rotamer_entry> const_rotamer_entry_shared_ptr;
typedef std::vector<rotamer_entry_shared_ptr> rotamer_entry_shared_ptr_vector;
typedef std::vector<const_rotamer_entry_shared_ptr> const_rotamer_entry_shared_ptr_vector;

inline const int get_combined_r(const int r1,const int r2,const int r3,const int r4){
	return r1
		+ (10 * r2)
		+ (100 * r3)
		+ (1000 * r4);
}

inline const int get_combined_pp_bin(const int phi, const int psi){
	return phi
		+ (100 * psi);
}

inline const int get_dih_deg_bin(const double ang_degrees, const double bin_size){
	//const double bin_tol = 0.00001;
	const double half_bin = bin_size / 2.0;
	const double max_ang = 360.0;
	const int max_bin = static_cast<int>(max_ang / bin_size) - 1;
	const double mod_ang = (ang_degrees + 180.0);
	int bin = static_cast<int>( (mod_ang + half_bin) /bin_size);
	//if (mod_ang >= 360 && mod_ang < 360.0 + bin_tol) mod_ang = 0.0;
	if (bin == max_bin+1) bin = 0;
	return bin;
}

inline const int get_dih_rad_bin(const double ang_radians, const double bin_size){
	return get_dih_deg_bin(UTILS::radians_to_degrees(ang_radians), UTILS::radians_to_degrees(bin_size));
}

class rotamer_entry {


	friend rotamer_entry_shared_ptr new_rotamer_entry(const POSE::residue_type type,
			const double phi,
			const double psi,
			const int r1,
			const int r2,
			const int r3,
			const int r4,
			const double prob,
			const double chi1,
			const double chi2,
			const double chi3,
			const double chi4,
			const double chi1_sig,
			const double chi2_sig,
			const double chi3_sig,
			const double chi4_sig,
			const double pp_bin_size_deg,
			const bool is_bb_dep);

private:

	rotamer_entry(const POSE::residue_type type,
			const double phi,
			const double psi,
			const int r1,
			const int r2,
			const int r3,
			const int r4,
			const double prob,
			const double chi1,
			const double chi2,
			const double chi3,
			const double chi4,
			const double chi1_sig,
			const double chi2_sig,
			const double chi3_sig,
			const double chi4_sig,
			const double pp_bin_size_deg,
			const bool is_bb_dep);
	//rotamer_entry();
	rotamer_entry(const rotamer_entry&);



public:

	const POSE::residue_type type;
	const double phi;
	const double psi;
	const int r1;
	const int r2;
	const int r3;
	const int r4;
	const double prob;
	const double chi1;
	const double chi2;
	const double chi3;
	const double chi4;
	const double chi1_sig;
	const double chi2_sig;
	const double chi3_sig;
	const double chi4_sig;
	const bool is_bb_dep;

	const double chi1_rad;
	const double chi2_rad;
	const double chi3_rad;
	const double chi4_rad;

	const double chi1_rad_sig;
	const double chi2_rad_sig;
	const double chi3_rad_sig;
	const double chi4_rad_sig;

	const int pp_bin_size_deg;

	const int phi_bin;
	const int psi_bin;

	const int combined_r;
	const int combined_pp_bin;



};

rotamer_entry_shared_ptr new_rotamer_entry(const POSE::residue_type type,
			const double phi,
			const double psi,
			const int r1,
			const int r2,
			const int r3,
			const int r4,
			const double prob,
			const double chi1,
			const double chi2,
			const double chi3,
			const double chi4,
			const double chi1_sig,
			const double chi2_sig,
			const double chi3_sig,
			const double chi4_sig,
			const double pp_bin_size_deg,
			const bool is_bb_dep);


}
}

#endif /* ROTAMER_ENTRY_H_ */
