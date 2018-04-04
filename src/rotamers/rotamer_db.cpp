/*
 * rotamer_db_interface.cpp
 *
 *  Created on: 29 Aug 2011
 *      Author: jmacdona
 */



#include "rotamer_db.h"

using namespace boost;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace std;

namespace PRODART {
namespace ROTAMERS {


rotamer_db::rotamer_db(){
	this->clear();
}


rotamer_db_shared_ptr new_rotamer_db(){

	boost::shared_ptr<rotamer_db> nat(new rotamer_db);
	return nat;
}


bool rotamer_db::load_gz_data( std::istream& input ){
	using namespace boost::iostreams;

	boost::iostreams::filtering_istreambuf in;



	in.push(gzip_decompressor());
	in.push(input);
	std::istream instr(&in);

	return this->load_data(instr);

}

void rotamer_db::clear(){
	all_rotamers.clear();
	all_rotamers.reserve(1000000);
	rt_pp_bin_r_bin_rotomers.clear();
	rt_pp_bin_rotomers_vec.clear();
	pp_bin_size_deg = 10.0;
	pp_bin_size_rad = UTILS::degrees_to_radians(pp_bin_size_deg);
	pp_ignore_val = 180;
	use_ignore_val = true;
	is_bb_dep = true;
}


void rotamer_db::index_rotamers(){

	rotamer_entry_shared_ptr_vector::iterator it;

	for (it = all_rotamers.begin(); it != all_rotamers.end(); it++){
		const residue_type rt = (*it)->type;
		const int pp_bin = (*it)->combined_pp_bin;
		const int r_bin = (*it)->combined_r;
		const_rotamer_entry_shared_ptr ptr = *it;
		rt_pp_bin_r_bin_rotomers[rt][pp_bin][r_bin] = ptr;
		rt_pp_bin_rotomers_vec[rt][pp_bin].push_back(ptr);
	}

}


const_rotamer_entry_shared_ptr rotamer_db::get_rotamer_entry(const POSE::residue_type& type,
			const int combined_pp_bin,
			const int combined_r_bin) const{
	const POSE::residue_type& rt (type.get_label3());

	rt_int_int_const_rotamer_entry_shared_ptr_unord_map::const_iterator it1 = rt_pp_bin_r_bin_rotomers.find(rt);
	if (it1 != rt_pp_bin_r_bin_rotomers.end()){
		int_int_const_rotamer_entry_shared_ptr_unord_map::const_iterator it2 = it1->second.find(combined_pp_bin);
		if (it2 != it1->second.end()){
			int_const_rotamer_entry_shared_ptr_unord_map::const_iterator it3 = it2->second.find(combined_r_bin);
			if (it3 != it2->second.end()){
				return it3->second;
			}

		}
	}

	return const_rotamer_entry_shared_ptr();
}

const_rotamer_entry_shared_ptr rotamer_db::get_rotamer_entry(const POSE::residue_type& type,
			const double phi,
			const double psi,
			const double chi1,
			const double chi2,
			const double chi3,
			const double chi4) const{
	//const POSE::residue_type& rt (type.get_label3());

	const int phi_bin = get_dih_rad_bin(phi, pp_bin_size_rad);
	const int psi_bin = get_dih_rad_bin(psi, pp_bin_size_rad);
	const int r1 = get_chi_r(type, 1, chi1);
	const int r2 = get_chi_r(type, 2, chi2);
	const int r3 = get_chi_r(type, 3, chi3);
	const int r4 = get_chi_r(type, 4, chi4);

	const int combined_pp_bin = get_combined_pp_bin(phi_bin, psi_bin);
	const int combined_r_bin = get_combined_r(r1, r2, r3, r4);

	return get_rotamer_entry(type, combined_pp_bin, combined_r_bin);

	/*
	rt_int_int_const_rotamer_entry_shared_ptr_unord_map::const_iterator it1 = rt_pp_bin_r_bin_rotomers.find(rt);
	if (it1 != rt_pp_bin_r_bin_rotomers.end()){
		int_int_const_rotamer_entry_shared_ptr_unord_map::const_iterator it2 = it1->second.find(combined_pp_bin);
		if (it2 != it1->second.end()){
			int_const_rotamer_entry_shared_ptr_unord_map::const_iterator it3 = it2->second.find(combined_r_bin);
			if (it3 != it2->second.end()){
				return it3->second;
			}

		}
	}

	return const_rotamer_entry_shared_ptr();
	*/
}


const const_rotamer_entry_shared_ptr_vector rotamer_db::get_rotamer_entry_vector(const POSE::residue_type& type,
			const int combined_pp_bin) const{
	const POSE::residue_type& rt (type.get_label3());
	rt_int_const_rotamer_entry_shared_ptr_vector_unord_map::const_iterator it1 = rt_pp_bin_rotomers_vec.find(rt);
	if (it1 != rt_pp_bin_rotomers_vec.end()){
		int_const_rotamer_entry_shared_ptr_vector_unord_map::const_iterator it2 = it1->second.find(combined_pp_bin);
		if (it2 != it1->second.end()){
			return it2->second;
		}
	}


	return const_rotamer_entry_shared_ptr_vector(0);
}

const const_rotamer_entry_shared_ptr_vector rotamer_db::get_rotamer_entry_vector(const POSE::residue_type& type,
			const double phi,
			const double psi) const{
	const int phi_bin = get_dih_rad_bin(phi, pp_bin_size_rad);
	const int psi_bin = get_dih_rad_bin(psi, pp_bin_size_rad);
	const int combined_pp_bin = get_combined_pp_bin(phi_bin, psi_bin);
	return get_rotamer_entry_vector(type, combined_pp_bin);
}

bool rotamer_db::load_data( std::istream& input ) {

	this->clear();

	string lineStr;


	long length, lineNum = 0 ;

	long not_parsed = 0;

	std::vector<std::string> SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		trim(lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t "), boost::token_compress_on );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 17 ){

				residue_type type(lexical_cast<string>(SplitVec[0]));
				const double phi = lexical_cast<double>(SplitVec[1]);
				const double psi = lexical_cast<double>(SplitVec[2]);
				//const long count = lexical_cast<double>(SplitVec[3]);
				const int r1 = lexical_cast<int>(SplitVec[4]);
				const int r2 = lexical_cast<int>(SplitVec[5]);
				const int r3 = lexical_cast<int>(SplitVec[6]);
				const int r4 = lexical_cast<int>(SplitVec[7]);
				const double prob = lexical_cast<double>(SplitVec[8]);
				const double chi1 = lexical_cast<double>(SplitVec[9]);
				const double chi2 = lexical_cast<double>(SplitVec[10]);
				const double chi3 = lexical_cast<double>(SplitVec[11]);
				const double chi4 = lexical_cast<double>(SplitVec[12]);
				const double chi1_sig = lexical_cast<double>(SplitVec[13]);
				const double chi2_sig = lexical_cast<double>(SplitVec[14]);
				const double chi3_sig = lexical_cast<double>(SplitVec[15]);
				const double chi4_sig = lexical_cast<double>(SplitVec[16]);

				if ( (phi != pp_ignore_val && psi != pp_ignore_val) || use_ignore_val == false){

					rotamer_entry_shared_ptr entry = new_rotamer_entry(type,
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
							is_bb_dep);

					all_rotamers.push_back(entry);
				}

				//cout << SplitVec[0] << " : " << SplitVec[16]  << endl;


			}
			else {
				not_parsed++;
			}
		}



	}

	this->index_rotamers();

	cout << "rotamer_db: loaded " << all_rotamers.size() << " rotamers" << endl;

	if (not_parsed != 0){
		cerr << "rotamer_db: WARNING: "  << not_parsed <<  " lines were not parsed due to format issues" << endl;
	}

	return true;
}


}
}

