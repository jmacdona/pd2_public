/*
 * ca_hbond.cpp
 *
 *  Created on: 1 Sep 2010
 *      Author: jmacdona
 */
#include "ca_hbond.h"
using namespace boost;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;
using namespace std;

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{



potential_shared_ptr new_ca_hbond(){
	potential_shared_ptr ptr(new ca_hbond());
	return ptr;
}

const int num_types = 4;
const int num_bins = 20;
const double dist_cutoff = 4.5;
const double min_dist = 3.0;
const double chi_min = (PI / 2.0);



ca_hbond::ca_hbond() {

	name_vector.clear();
	name_vector.push_back(potentials_name("ca_hb_other"));
	name_vector.push_back(potentials_name("ca_hb_alpha"));
	name_vector.push_back(potentials_name("ca_hb_sr_beta"));
	name_vector.push_back(potentials_name("ca_hb_lr_beta"));

	hbond_energy.resize(num_types, 0);
	N_O_dist_bin_score.resize(num_types);
	chi_i_bin_score.resize(num_types);
	chi_j_bin_score.resize(num_types);
	tau_bin_score.resize(num_types);
	for (int i = 0; i < num_types; i++) {
		N_O_dist_bin_score[i].resize(num_bins, 0);
		chi_i_bin_score[i].resize(num_bins, 0);
		chi_j_bin_score[i].resize(num_bins, 0);
		tau_bin_score[i].resize(num_bins, 0);
	}


}

int ca_hbond::getHbond_type(PRODART::POSE::four_state_sec_struct phipsiombin_i,
		PRODART::POSE::four_state_sec_struct phipsiombin_j ,
		int seq_sep) const{

	if (seq_sep < 3){
		return 4; //no hbond
	}
	else if (phipsiombin_i == PRODART::POSE::ss4_CIS
			|| phipsiombin_j == PRODART::POSE::ss4_CIS){
		return 4; //no hbond
	}
	else if ((phipsiombin_i == PRODART::POSE::ss4_HELIX
			|| phipsiombin_j == PRODART::POSE::ss4_HELIX) && seq_sep == 3){ //seq_sep ==3 due to pseudo-atom numbering
		return 1; //alpha helix
	}
	else if ((phipsiombin_i == PRODART::POSE::ss4_STRAND
			|| phipsiombin_j == PRODART::POSE::ss4_STRAND)
			&& seq_sep <= 5
			&& seq_sep > 2){

		return 2; // short range beta sheet

	}
	else if ((phipsiombin_i == PRODART::POSE::ss4_STRAND
			&& phipsiombin_j == PRODART::POSE::ss4_STRAND) && seq_sep > 5){

		return 3; // long range beta sheet

	}
	else {
		return 0;
	}

}


bool ca_hbond::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:ca_hbond");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	return true;
}

/*
bool ca_hbond::provides(const potentials_name& query_name) const{

	potentials_name_vector::const_iterator iter;
	for (iter = this->name_vector.begin(); iter != this->name_vector.end(); iter++){
		if (query_name == *iter){
			return true;
		}
	}

	return false;
}
*/

std::istream& ca_hbond::load_data( std::istream& input ){


	string lineStr;


	long length, lineNum = 0 ;


	string_vector SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 8 ){

				string paraName = SplitVec[1];
				trim(paraName);

				const int bin = lexical_cast<int>(SplitVec[2]);

				if ( paraName.compare("N_O_dist") == 0 ) {
					/*
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double par_beta = lexical_cast<double>(SplitVec[6]);
					const double long_beta = lexical_cast<double>(SplitVec[7]);
					other_N_O_dist_bin_score[bin] = other;
					alpha_N_O_dist_bin_score[bin] = alpha;
					short_beta_N_O_dist_bin_score[bin] = par_beta;
					long_beta_N_O_dist_bin_score[bin] = long_beta;
					*/

					for (int i = 0; i < num_types; i++){
						N_O_dist_bin_score[i][bin] = lexical_cast<double>(SplitVec[4+i]);
					}




				}
				else if ( paraName.compare("chi_i") == 0 ){
					/*
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double par_beta = lexical_cast<double>(SplitVec[6]);
					const double long_beta = lexical_cast<double>(SplitVec[7]);

					other_chi_i_bin_score[bin] = other;
					alpha_chi_i_bin_score[bin] = alpha;
					short_beta_chi_i_bin_score[bin] = par_beta;
					long_beta_chi_i_bin_score[bin] = long_beta;
					*/

					for (int i = 0; i < num_types; i++){
						chi_i_bin_score[i][bin] = lexical_cast<double>(SplitVec[4+i]);
					}

				}
				else if ( paraName.compare("chi_j") == 0 ){
					/*
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double par_beta = lexical_cast<double>(SplitVec[6]);
					const double long_beta = lexical_cast<double>(SplitVec[7]);

					other_chi_j_bin_score[bin] = other;
					alpha_chi_j_bin_score[bin] = alpha;
					short_beta_chi_j_bin_score[bin] = par_beta;
					long_beta_chi_j_bin_score[bin] = long_beta;
					*/

					for (int i = 0; i < num_types; i++){
						chi_j_bin_score[i][bin] = lexical_cast<double>(SplitVec[4+i]);
					}


				}
				else if ( paraName.compare("tau") == 0 ){
					/*
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double par_beta = lexical_cast<double>(SplitVec[6]);
					const double long_beta = lexical_cast<double>(SplitVec[7]);

					other_tau_bin_score[bin] = other;
					alpha_tau_bin_score[bin] = alpha;
					short_beta_tau_bin_score[bin] = par_beta;
					long_beta_tau_bin_score[bin] = long_beta;
					*/

					for (int i = 0; i < num_types; i++){
						tau_bin_score[i][bin] = lexical_cast<double>(SplitVec[4+i]);
					}


				}
				else if ( paraName.compare("hbond_energy") == 0 ){
					/*
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double par_beta = lexical_cast<double>(SplitVec[6]);
					const double long_beta = lexical_cast<double>(SplitVec[7]);

					other_hb_energy = other;
					alpha_hb_energy = alpha;
					short_beta_hb_energy = par_beta;
					long_beta_hb_energy = long_beta;
					*/

					for (int i = 0; i < num_types; i++){
						hbond_energy[i] = lexical_cast<double>(SplitVec[4+i]);
					}

					/*
					cout << "hbond_energy\t"
						 <<  other_hb_energy << "\t"
						 << alpha_hb_energy << "\t"
						 << beta_hb_energy << "\t"
						 <<	endl;
					*/


				}
				else {
                    cout << "ERROR - unknown parameter name: " << paraName << endl;
                    //fatalError = true;
                    //cout << paraValue << endl;
                }


			}
		}
	}




	return input;

}


inline int ca_hbond::get_N_O_dist_bin(double dist) const{

	const double inc = 1.5 / 20.0;
	return static_cast<int>((dist - min_dist)/inc);

}

inline int ca_hbond::get_chi_i_bin(double chi_i) const{
	const double inc = (PI / 2.0) / 20.0;
	return static_cast<int>((chi_i - chi_min )/inc);
}

inline int ca_hbond::get_chi_j_bin(double chi_j) const{
	const double inc = (PI / 2.0) / 20.0;
	return static_cast<int>((chi_j - chi_min )/inc);
}

inline int ca_hbond::get_tau_bin(double tau) const{
	const double inc = (PI * 2.0) / 20.0;
	return static_cast<int>((tau + PI )/inc);
}


double ca_hbond::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	PRODART::POSE::META::nb_ele_vector& pair_list = ca_meta_dat->get_NO_pair_list();
	PRODART::POSE::four_state_sec_struct_vector& conf_class = ca_meta_dat->get_conf_class();
	double_vector enrg_components(4 , 0);

	PRODART::POSE::META::nb_ele_vector::iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		//cout << "\ndb: 1 ";
		if (iter->dist < dist_cutoff){
			//cout << "db: 2 ";
			const int seq_sep = iter->seq_sep;
			PRODART::POSE::four_state_sec_struct ss1 = conf_class[iter->res_num1];
			PRODART::POSE::four_state_sec_struct ss2 = conf_class[iter->res_num2];
			//cout << "db: " << iter->atype1.get_label() << " " << iter->atype2.get_label() << " ";
			if ((iter->atype1 == ca_pose_meta::pseudo_N() &&  iter->atype2 == ca_pose_meta::pseudo_O())
					|| (iter->atype1 == ca_pose_meta::pseudo_O() &&  iter->atype2 == ca_pose_meta::pseudo_N()) ){
				//cout << "db: 3 ";
				const int hb_type = this->getHbond_type(ss1, ss2, seq_sep);
				ca_hbond_collection_shared_ptr n_collect;
				ca_hbond_collection_shared_ptr o_collect;

				if ((iter->atype1 == ca_pose_meta::pseudo_N()
						&&  iter->atype2 == ca_pose_meta::pseudo_O())){
					//cout << "db: a\n";
					n_collect = ca_meta_dat->get_hbond_collection(iter->atom1_ptr);
					o_collect = ca_meta_dat->get_hbond_collection(iter->atom2_ptr);
				}
				else {
					//cout << "db: b\n";
					n_collect = ca_meta_dat->get_hbond_collection(iter->atom2_ptr);
					o_collect = ca_meta_dat->get_hbond_collection(iter->atom1_ptr);
				}


				if (hb_type < num_types){
					//cout << "db: 4 ";
					//for some odd reason chi_i and chi_j are swapped in the param file
					const double chi_j = angle(o_collect->N_i->get_coords(), o_collect->O_i->get_coords(), n_collect->N_i->get_coords());
					const double chi_i = angle(n_collect->O_i->get_coords(), n_collect->N_i->get_coords(), o_collect->O_i->get_coords());
					if (chi_i > chi_min && chi_j > chi_min ){
						//cout << "db: 5 ";
						const double tau = dihedral(o_collect->CA_i->get_coords(),
								o_collect->N_i->get_coords(),
								n_collect->N_i->get_coords(),
								n_collect->CA_i->get_coords());

						const double this_dist = iter->dist < min_dist ? min_dist : iter->dist;

						const int N_O_dist_bin = get_N_O_dist_bin(this_dist);
						const int chi_i_bin = get_chi_i_bin(chi_i);
						const int chi_j_bin = get_chi_j_bin(chi_j);
						const int tau_bin = get_tau_bin(tau);

						double this_energy = 0;

						this_energy += N_O_dist_bin_score[hb_type][N_O_dist_bin];
						this_energy += chi_i_bin_score[hb_type][chi_i_bin];
						this_energy += chi_j_bin_score[hb_type][chi_j_bin];
						this_energy += tau_bin_score[hb_type][tau_bin];
						this_energy += hbond_energy[hb_type];

						enrg_components[hb_type] += this_energy;

						/*
						cout << "dist " << N_O_dist_bin << " " << this_dist << " " << std::flush;
						cout << "chi_i " << chi_i_bin << " " << chi_i << " " << std::flush;
						cout << "chi_j " << chi_j_bin << " " << chi_j << " " << std::flush;
						cout << "tau " << tau_bin << " " << tau << " " << std::flush;
						cout << "NO_enrg " <<  N_O_dist_bin_score[hb_type][N_O_dist_bin] << " ";
						cout << "chi_i_enrg " <<  chi_i_bin_score[hb_type][chi_i_bin] << " ";
						cout << "chi_j_enrg " <<  chi_j_bin_score[hb_type][chi_j_bin] << " ";
						cout << "tau_enrg " <<  tau_bin_score[hb_type][tau_bin] << " ";
						cout << iter->res_num1 <<  " " << iter->res_num2 << " " << std::flush;
						cout << "seq_sep: " << iter->seq_sep << " " << std::flush;
						cout << "hbond_type: " << hb_type << " " << std::flush;
						cout << "this_energy: " << this_energy << " " << std::flush;
						cout << iter->atom1_ptr->get_residue()->get_pdb_residue_index() <<  " "
								<< iter->atom2_ptr->get_residue()->get_pdb_residue_index() << " "
								<< std::endl;
						 */

					}
				}


			}


		}
	}

	double total_score = 0;
	for (int i = 0; i < num_types; i++) {
		total_score += energies_map.add_energy_component(name_vector[i], enrg_components[i]);
	}
	return total_score;
}

double ca_hbond::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	return this->get_energy(pose_meta_, energies_map);
}


}
}
}
}
