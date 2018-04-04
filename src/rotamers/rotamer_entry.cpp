/*
 * rotamer_entry.cpp
 *
 *  Created on: 28 Aug 2011
 *      Author: jmacdona
 */
#include "rotamer_entry.h"

using namespace std;


namespace PRODART {
namespace ROTAMERS {

/*
rotamer_entry::rotamer_entry(){

}

rotamer_entry::rotamer_entry(const rotamer_entry&){

}
*/

/*
const int get_dih_deg_bin(const double ang, const double bin_size){
	//const double bin_tol = 0.00001;
	const double half_bin = bin_size / 2.0;
	const double max_ang = 360.0;
	const int max_bin = static_cast<int>(max_ang / bin_size) - 1;
	const double mod_ang = (ang + 180.0);
	int bin = static_cast<int>( (mod_ang + half_bin) /bin_size);
	//if (mod_ang >= 360 && mod_ang < 360.0 + bin_tol) mod_ang = 0.0;
	if (bin == max_bin+1) bin = 0;
	return bin;
}
*/

rotamer_entry::rotamer_entry(const POSE::residue_type type_,
		const double phi_,
		const double psi_,
		const int r1_,
		const int r2_,
		const int r3_,
		const int r4_,
		const double prob_,
		const double chi1_,
		const double chi2_,
		const double chi3_,
		const double chi4_,
		const double chi1_sig_,
		const double chi2_sig_,
		const double chi3_sig_,
		const double chi4_sig_,
		const double pp_bin_size_deg_,
		const bool is_bb_dep_) : type(type_),
		phi(phi_),
		psi(psi_),
		r1(r1_),
		r2(r2_),
		r3(r3_),
		r4(r4_),
		prob(prob_),
		chi1(chi1_),
		chi2(chi2_),
		chi3(chi3_),
		chi4(chi4_),
		chi1_sig(chi1_sig_),
		chi2_sig(chi2_sig_),
		chi3_sig(chi3_sig_),
		chi4_sig(chi4_sig_),
		is_bb_dep(is_bb_dep_),
		chi1_rad(UTILS::degrees_to_radians(chi1_)),
		chi2_rad(UTILS::degrees_to_radians(chi2_)),
		chi3_rad(UTILS::degrees_to_radians(chi3_)),
		chi4_rad(UTILS::degrees_to_radians(chi4_)),
		chi1_rad_sig(UTILS::degrees_to_radians(chi1_sig_)),
		chi2_rad_sig(UTILS::degrees_to_radians(chi2_sig_)),
		chi3_rad_sig(UTILS::degrees_to_radians(chi3_sig_)),
		chi4_rad_sig(UTILS::degrees_to_radians(chi4_sig_)),
		pp_bin_size_deg(pp_bin_size_deg_),
		phi_bin(get_dih_deg_bin(phi_, pp_bin_size_deg_)),
		psi_bin(get_dih_deg_bin(psi_, pp_bin_size_deg_)),
		combined_r(get_combined_r(r1_, r2_,r3_,r4_)),
		combined_pp_bin(get_combined_pp_bin(phi_bin, psi_bin)){

	//cout << phi << " " << phi_bin << endl;

}

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
			const bool is_bb_dep){
	boost::shared_ptr<rotamer_entry> re(new rotamer_entry(type,
			phi,
			psi,
			r1,
			r2,
			r3,
			r4,
			prob,
			chi1,
			chi2,
			chi3,
			chi4,
			chi1_sig,
			chi2_sig,
			chi3_sig,
			chi4_sig,
			pp_bin_size_deg,
			is_bb_dep));
	return re;
}

















}
}
