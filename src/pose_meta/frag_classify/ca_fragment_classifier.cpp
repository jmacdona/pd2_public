/*
 * ca_fragment_classifier.cpp
 *
 *  Created on: 25 Aug 2010
 *      Author: jmacdona
 */
#include "ca_fragment_classifier.h"

using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;

using namespace std;

using boost::split;
using boost::is_any_of;
using boost::trim;
using boost::lexical_cast;

namespace PRODART {
namespace POSE {
namespace META {
namespace FRAG_CLASS {

typedef std::vector< std::string >  split_vector_type;


t_frag_count_element::t_frag_count_element(void){
	//cout << "t_IBBB_count_element constructor" << endl;
	ss_bin_count.clear();
	frag_number_count.clear();
	total_count = 0;
}

t_frag_count_element::~t_frag_count_element(void){

}

namespace {
fragment_classifier_shared_ptr instance = fragment_classifier_shared_ptr();
}

fragment_classifier_shared_ptr ca_fragment_classifier::MainInstance(){
	if (!instance){
		//fragment_classifier_shared_ptr ptr(new ca_fragment_classifier());
		//instance = ptr;//new ca_fragment_classifier(); //ptr;//fragment_classifier_shared_ptr(new ca_fragment_classifier());
		instance = new_ca_fragment_classifier();
	}
	return instance;
}



fragment_classifier_shared_ptr new_ca_fragment_classifier(){
	fragment_classifier_shared_ptr ptr(new ca_fragment_classifier());
	return ptr;
}

const double ca_fragment_classifier::four_CA_bin_size = 0.3;

ca_fragment_classifier::ca_fragment_classifier(){
	const string db_path = PRODART::ENV::get_option_value<string>("database:path:ca_fragment_classifier_db");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->loadData(inputdb);
	inputdb.close();
}

bool ca_fragment_classifier::classify_fragments (const pose_shared_ptr _pose,
		const PRODART::POSE::META::ca_pose_meta_shared_ptr _pose_meta) const{

	frag4_vector &fragments = _pose_meta->get_fragments();
	vector3d CA_vecs[4];

	for (unsigned int i = 0 ; i < fragments.size(); i++){
		frag4_element &ele = fragments[i];
		for (int v = 0 ; v < 4; v++){
			CA_vecs[v] = ele.CA_atoms[v]->get_coords();
		}



		const double CA_tau = dihedral(CA_vecs[0], CA_vecs[1], CA_vecs[2], CA_vecs[3]);
		const double CA_omega1 = angle(CA_vecs[0], CA_vecs[1], CA_vecs[2]);
		const double CA_omega2 = angle(CA_vecs[1], CA_vecs[2], CA_vecs[3]);

		ele.tau = CA_tau;
		ele.omega1 = CA_omega1;
		ele.omega2 = CA_omega2;

		const vector3d v_i_minus1 = CA_vecs[1] - CA_vecs[0];
		const vector3d v_i = CA_vecs[2] - CA_vecs[1];
		const vector3d v_i_plus1 = CA_vecs[3] - CA_vecs[2];
		//const vector3d v_p = CA_vecs[3] - CA_vecs[1];

		double chi = (v_i_minus1 * v_i).dot(v_i_plus1);

		if (chi == 0){
			chi = 1;
		}
		else {
			chi = chi / fabs(chi);
		}

		const double d_m1_p1 = (CA_vecs[0] - CA_vecs[2]).mod();//meta_dat.getCADist(overall_res_i_index - 1, overall_res_i_index + 1); //(CA_vecs[0] - CA_vecs[2]).mod();
		const double d_0_p2 = (CA_vecs[1] - CA_vecs[3]).mod();//meta_dat.getCADist(overall_res_i_index, overall_res_i_index + 2); //(CA_vecs[1] - CA_vecs[3]).mod();
		const double d_m1_p2 = chi * (CA_vecs[0] - CA_vecs[3]).mod(); //chi * meta_dat.getCADist(overall_res_i_index - 1, overall_res_i_index + 2); //(CA_vecs[0] - CA_vecs[3]).mod();

		const int d_m1_p1_bin = static_cast<int>( d_m1_p1 / four_CA_bin_size );
		const int d_0_p2_bin = static_cast<int>( d_0_p2 / four_CA_bin_size );
		const int d_m1_p2_bin = static_cast<int>( d_m1_p2 / four_CA_bin_size );

		const int overall_key = static_cast<int>(chi) * (10000 * abs(d_m1_p2_bin)
						+ 100 * d_0_p2_bin
						+  d_m1_p1_bin);

		ele.IBBB_key = overall_key;

		//find approximate closest key using lower_bound to save time
		int_int_map::const_iterator ibbb_iter = IBBB_to_frag_number_map.lower_bound(overall_key);
		int frag_number = 1;
		if (ibbb_iter ==  IBBB_to_frag_number_map.end()){
			const int closest_key = this->findClosestKey(overall_key);
			int_int_map::const_iterator ibbb_iter2 = IBBB_to_frag_number_map.find(closest_key);
			if (ibbb_iter2 !=  IBBB_to_frag_number_map.end()){
				frag_number = ibbb_iter2->second;
			}
			else {
				cerr << "\nWARNING: ca_fragment_classifier: IBBB fragment assignment failure!\n";
				cerr << "CA[1] residue index:\t" << ele.CA_1_residue_number << "\n";
				cerr << "overall_key:\t" << overall_key << "\n";
			}

			//return false;
		}
		else {
			frag_number = ibbb_iter->second; //IBBB_to_frag_number_map[overall_key];
		}

		int_int_umap::const_iterator sec_iter = frag_number_to_SStype_map.find(frag_number);
		int ss_type = 0;
		if (sec_iter == frag_number_to_SStype_map.end()){
			cerr << "\nWARNING: ca_fragment_classifier: sec struct fragment assignment failure!\n" << "\n";
			//return false;
		}
		else {
			ss_type = sec_iter->second; //frag_number_to_SStype_map[frag_number];
		}

		ele.frag_type_num = frag_number;
		ele.frag_ss_class = static_cast<PRODART::POSE::four_state_sec_struct>(ss_type);

		_pose_meta->set_residue_conf_class(ele.CA_1_residue_number, ele.frag_ss_class);
		_pose_meta->set_frag_type_num(ele.CA_1_residue_number, ele.frag_type_num);
		if (ele.frag_pos == fragNTERM || ele.frag_pos == fragNandCTERM){
			_pose_meta->set_residue_conf_class(ele.CA_1_residue_number - 1, ele.frag_ss_class);
			_pose_meta->set_frag_type_num(ele.CA_1_residue_number - 1, ele.frag_type_num);
		}
		if (ele.frag_pos == fragCTERM || ele.frag_pos == fragNandCTERM){
			_pose_meta->set_residue_conf_class(ele.CA_1_residue_number + 1, ele.frag_ss_class);
			_pose_meta->set_residue_conf_class(ele.CA_1_residue_number + 2, ele.frag_ss_class);
			_pose_meta->set_frag_type_num(ele.CA_1_residue_number + 1, ele.frag_type_num);
			_pose_meta->set_frag_type_num(ele.CA_1_residue_number + 2, ele.frag_type_num);
		}

		const vector3d N_prime_coords = (CA_vecs[1] + CA_vecs[2]) / 2;
		const vector3d o_prime_unnorm = v_i * v_i_plus1;
		const vector3d o_prime_norm = o_prime_unnorm / o_prime_unnorm.mod();
		const vector3d o_prime_coords = N_prime_coords + o_prime_norm;

		_pose->set_atom_coords(N_prime_coords,ca_pose_meta::pseudo_N(), ele.CA_1_residue_number);
		_pose->set_atom_coords(o_prime_coords,ca_pose_meta::pseudo_O(), ele.CA_1_residue_number);

		if (ele.frag_pos == fragNTERM || ele.frag_pos == fragNandCTERM){
			const vector3d N_prime_coords0 = (CA_vecs[0] + CA_vecs[1]) / 2;
			const vector3d o_prime_unnorm0 = (CA_vecs[1] - CA_vecs[0]) * (CA_vecs[2] - CA_vecs[1]);
			const vector3d o_prime_norm0 = o_prime_unnorm0 / o_prime_unnorm0.mod();
			const vector3d o_prime_coords0 = N_prime_coords0 + o_prime_norm0;

			_pose->set_atom_coords(N_prime_coords0, ca_pose_meta::pseudo_N(), ele.CA_1_residue_number-1);
			_pose->set_atom_coords(o_prime_coords0, ca_pose_meta::pseudo_O(), ele.CA_1_residue_number-1);

		}


	}
	return true;
}

inline double ca_fragment_classifier::IBBB_key_distance(const int key1, const int key2) const{

	const int abs_key1 = abs(key1);
	const int abs_key2 = abs(key2);

	const int sign_key1 = key1 / abs_key1;
	const int sign_key2 = key2 / abs_key2;

	/*
	const int overall_key = static_cast<int>(chi) * (10000 * abs(d_m1_p2_bin)
							+ 100 * d_0_p2_bin
							+  d_m1_p1_bin);
	*/
	const int key1_d_m1_p2_bin = sign_key1 * (abs_key1 / 10000);
	const int key1_d_0_p2_bin = (abs_key1 - abs(key1_d_m1_p2_bin * 10000))
								/ 100;
	const int key1_d_m1_p1_bin = (abs_key1 - abs(key1_d_m1_p2_bin * 10000)
								- (key1_d_0_p2_bin * 100));

	const int key2_d_m1_p2_bin = sign_key2 * (abs_key2 / 10000);
	const int key2_d_0_p2_bin = (abs_key2 - abs(key2_d_m1_p2_bin * 10000))
								/ 100;
	const int key2_d_m1_p1_bin = (abs_key2 - abs(key2_d_m1_p2_bin * 10000)
								- (key2_d_0_p2_bin * 100));

	/*
	cout << key1 << "\t"
		 << key1_d_m1_p2_bin << "\t"
		 << key1_d_0_p2_bin << "\t"
		 << key1_d_m1_p1_bin << "\t"
		 << endl;

	cout << key2 << "\t"
		 << key2_d_m1_p2_bin << "\t"
		 << key2_d_0_p2_bin << "\t"
		 << key2_d_m1_p1_bin << "\t"
		 << endl;
	*/
	double distance = sqrt((key1_d_m1_p2_bin - key2_d_m1_p2_bin) * (key1_d_m1_p2_bin - key2_d_m1_p2_bin)
					+ (key1_d_0_p2_bin - key2_d_0_p2_bin) * (key1_d_0_p2_bin - key2_d_0_p2_bin)
					+ (key1_d_m1_p1_bin - key2_d_m1_p1_bin) * (key1_d_m1_p1_bin - key2_d_m1_p1_bin));

	/*
	cout << distance << endl;
	*/

	return distance;
}

inline int ca_fragment_classifier::findClosestKey(const int search_key) const{
	int_int_map::const_iterator iter;
	double best_dist = FLT_MAX;
	int best_key = 0;
	//int_int_umap::const_iterator lb = IBBB_to_frag_number_map
	int_int_map::const_iterator lb = IBBB_to_frag_number_map.lower_bound(search_key);
	if (lb != IBBB_to_frag_number_map.end()){
		return lb->first;
	}
	else {
		int_int_map::const_iterator ub = IBBB_to_frag_number_map.upper_bound(search_key);
		if (ub != IBBB_to_frag_number_map.end()){
			return ub->first;
		}
		else {

			if (search_key <= 0){
				int_int_map::const_iterator fl =  IBBB_to_frag_number_map.begin();
				if (fl != IBBB_to_frag_number_map.end()){
					return fl->first;
				}
			}
			else if (search_key >= 0){
				int_int_map::const_reverse_iterator fl =  IBBB_to_frag_number_map.rbegin();
				if (fl != IBBB_to_frag_number_map.rend()){
					return fl->first;
				}
			}

			cerr << "WARNING: ca_fragment_classifier::findClosestKey : reached last resort key: " << search_key << endl;


			//! change this to better search algorithm
			for ( iter = IBBB_to_frag_number_map.begin(); iter != IBBB_to_frag_number_map.end(); iter++ ){
				int key = iter->first;
				double dist = this->IBBB_key_distance(search_key, key);

				if (dist < best_dist ){
					best_dist = dist;
					best_key = key;
				}

			}
			return best_key;
		}
	}
}

void ca_fragment_classifier::fillGapsInIBBB_map(){

	if (IBBB_to_frag_number_map.size() == 0) return;

	int_int_map temp_map;
	temp_map.clear();

	//const int first_key = (IBBB_to_SStype_map.begin())->first;
	//const int last_key = (--(IBBB_to_SStype_map.end()))->first;

	//!TODO rebuild db with larger max bins
	const int max_d_m1_p1 = static_cast<int>(2.0*(4.0 / four_CA_bin_size));
	const int max_d_0_p2 = static_cast<int>(2.0*(4.0 / four_CA_bin_size));
	const int max_d_m1_p2 = static_cast<int>(3.0*(4.0 / four_CA_bin_size));
	const int min_bin = static_cast<int>((3.0 / four_CA_bin_size));

	for (int d_m1_p1_bin = min_bin; d_m1_p1_bin < max_d_m1_p1; d_m1_p1_bin++){
		for (int d_0_p2_bin = min_bin; d_0_p2_bin < max_d_0_p2; d_0_p2_bin++){
			for (int d_m1_p2_bin = min_bin; d_m1_p2_bin < max_d_m1_p2; d_m1_p2_bin++){
				/*
				const double d_m1_p1 = (CA_vecs[0] - CA_vecs[2]).mod();
				const double d_0_p2 = (CA_vecs[1] - CA_vecs[3]).mod();
				const double d_m1_p2 = chi * (CA_vecs[0] - CA_vecs[3]).mod();

				const int d_m1_p1_bin = static_cast<int>( d_m1_p1 / four_CA_bin_size );
				const int d_0_p2_bin = static_cast<int>( d_0_p2 / four_CA_bin_size );
				const int d_m1_p2_bin = static_cast<int>( d_m1_p2 / four_CA_bin_size );
				*/

				const int pos_key = 1 * (10000 * abs(d_m1_p2_bin)
										+ 100 * d_0_p2_bin
										+  d_m1_p1_bin);
				const int neg_key = -1 * (10000 * abs(d_m1_p2_bin)
										+ 100 * d_0_p2_bin
										+  d_m1_p1_bin);

				const int closest_pos_key = this->findClosestKey(pos_key);
				const int closest_neg_key = this->findClosestKey(neg_key);

				if (closest_pos_key != pos_key){
					const int ss_type = IBBB_to_frag_number_map[closest_pos_key];
					temp_map[pos_key] = ss_type;
				}
				if (closest_neg_key != neg_key){
					const int ss_type = IBBB_to_frag_number_map[closest_neg_key];
					temp_map[neg_key] = ss_type;
				}

			}
		}
	}

	int_int_map::iterator iter;

	for (iter = temp_map.begin(); iter != temp_map.end(); iter++){
		IBBB_to_frag_number_map[iter->first] = iter->second;
	}


}


istream& ca_fragment_classifier::loadData( istream& input){

	string lineStr;

	long length, lineNum = 0 ;

	IBBB_SS_count.clear();
	frag_number_to_SStype_map.clear();
	IBBB_to_frag_number_map.clear();

	split_vector_type SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 2 ){

				string paraName = SplitVec[0];
				trim(paraName);

				if ( paraName.compare("4FRAG_NUMBER_MAP") == 0 ){
					const int IBBB_bin = lexical_cast<int>(SplitVec[1]);
					const int frag_num = lexical_cast<int>(SplitVec[2]);
					IBBB_to_frag_number_map[IBBB_bin] = frag_num;
				}
				else if ( paraName.compare("4FRAG_SS_MAP") == 0 ){
					const int frag_num = lexical_cast<int>(SplitVec[1]);
					const int ss_type = lexical_cast<int>(SplitVec[2]);
					frag_number_to_SStype_map[frag_num] = ss_type;

				}
				else {
                    cout << "ERROR - unknown parameter name: " << paraName << endl;
                    //fatalError = true;
                    //cout << paraValue << endl;
                }


			}
		}
	}

	fillGapsInIBBB_map();
	return input;


}



}
}
}
}
