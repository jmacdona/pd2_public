/*
 * ca_density.cpp
 *
 *  Created on: 21 Oct 2013
 *      Author: jmacdona
 */


#include "ca_density.h"

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


potential_shared_ptr new_ca_density(){
	potential_shared_ptr ptr(new ca_density());
	return ptr;
}

ca_density::ca_density() : cutoff(12.0), min_seq_sep(10), max_density(99) {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_density"));

	alpha_counts.clear();
	beta_counts.clear();
	other_counts.clear();
	all_counts.clear();
	total_alpha = 0;
	total_beta = 0;
	total_other = 0;
	total_all = 0;

	bg_alpha_counts.clear();
	bg_beta_counts.clear();
	bg_other_counts.clear();
	bg_all_counts.clear();
	bg_total_alpha = 0;
	bg_total_beta = 0;
	bg_total_other = 0;
	bg_total_all = 0;

	alpha_pot.clear();
	beta_pot.clear();
	other_pot.clear();
	all_pot.clear();

}



int_vector ca_density::get_ca_counts(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_, const double cutoff_) const{
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	int_vector rtn_vec(pose_->get_residue_count(), 0);

	//int_vector alpha_rtn_vec(pose_->get_residue_count(), 0);
	//int_vector beta_rtn_vec(pose_->get_residue_count(), 0);
	//int_vector other_rtn_vec(pose_->get_residue_count(), 0);


	const nb_ele_vector& pair_list = ca_meta_dat->get_ca_pair_list();
	nb_ele_vector::const_iterator iter;

	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		if (iter->dist < cutoff_ && iter->seq_sep >= min_seq_sep){

			rtn_vec[iter->res_num1]++;
			rtn_vec[iter->res_num2]++;


			/*
			if (ca_meta_dat->get_residue_sec_struct(iter->res_num1) == ss3_HELIX){
				alpha_rtn_vec[iter->res_num1]++;
			}
			else if (ca_meta_dat->get_residue_sec_struct(iter->res_num1) == ss3_STRAND){
				beta_rtn_vec[iter->res_num1]++;
			}
			else {
				other_rtn_vec[iter->res_num1]++;
			}
			*/

			/*
			if (ca_meta_dat->get_residue_sec_struct(iter->res_num2) == ss3_HELIX){
				alpha_rtn_vec[iter->res_num1]++;
			}
			else if (ca_meta_dat->get_residue_sec_struct(iter->res_num2) == ss3_STRAND){
				beta_rtn_vec[iter->res_num2]++;
			}
			else {
				other_rtn_vec[iter->res_num2]++;
			}
			*/

		}
	}
	return rtn_vec;//boost::make_tuple(alpha_rtn_vec, beta_rtn_vec, other_rtn_vec);
}



double ca_density::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	int_vector counts = this->get_ca_counts(pose_meta_, cutoff);

	const const_pose_shared_ptr _pose = pose_meta_->get_pose();
	double total_score = 0;


	for (int i = 0 ; i < (int)counts.size(); i++){

		if (counts[i] <= max_density){
			if (pose_meta_->get_residue_sec_struct(i) == ss3_HELIX){
				total_score += alpha_pot.lower_bound(counts[i])->second;
			}
			else if (pose_meta_->get_residue_sec_struct(i) == ss3_STRAND){
				total_score += beta_pot.lower_bound(counts[i])->second;
			}
			else {
				total_score += other_pot.lower_bound(counts[i])->second;
			}
		}
		else {
			total_score += 100;
		}
	}

	return energies_map.add_energy_component(this->name_vector[0], total_score);
}

//TODO
double ca_density::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	return get_energy(_pose_meta, energies_map);
}


bool ca_density::init(){
	const string db_path = PRODART::ENV::get_option_value<string>("database:path:ca_density");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	return true;
}



bool ca_density::training_pose(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_){
	int_vector counts = this->get_ca_counts(pose_meta_, cutoff);
	for (int i = 0 ; i < (int)counts.size(); i++){
		all_counts[counts[i]] = all_counts[counts[i]] + 1;
		total_all += 1;
		if (pose_meta_->get_residue_sec_struct(i) == ss3_HELIX){
			alpha_counts[counts[i]] = alpha_counts[counts[i]] + 1;
			total_alpha += 1;
		}
		else if (pose_meta_->get_residue_sec_struct(i) == ss3_STRAND){
			beta_counts[counts[i]] = beta_counts[counts[i]] + 1;
			total_beta += 1;
		}
		else {
			other_counts[counts[i]] = other_counts[counts[i]] + 1;
			total_other += 1;
		}
	}

	return true;
}

bool ca_density::background_training_pose(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_){
	int_vector counts = this->get_ca_counts(pose_meta_, cutoff);
	for (int i = 0 ; i < (int)counts.size(); i++){
		bg_all_counts[counts[i]] = bg_all_counts[counts[i]] + 1;
		bg_total_all += 1;
		if (pose_meta_->get_residue_sec_struct(i) == ss3_HELIX){
			bg_alpha_counts[counts[i]] = bg_alpha_counts[counts[i]] + 1;
			bg_total_alpha += 1;
		}
		else if (pose_meta_->get_residue_sec_struct(i) == ss3_STRAND){
			bg_beta_counts[counts[i]] = bg_beta_counts[counts[i]] + 1;
			bg_total_beta += 1;
		}
		else {
			bg_other_counts[counts[i]] = bg_other_counts[counts[i]] + 1;
			bg_total_other += 1;
		}
	}

	return true;
}

bool ca_density::calc_potentials_using_training(){

	const double val = 1.0;
	//add pseudo counts
	for (int i = 0 ; i < 100; i++){
		alpha_counts[i] = alpha_counts[i] + val;
		total_alpha += val;

		beta_counts[i] = beta_counts[i] + val;
		total_beta += val;

		other_counts[i] = other_counts[i] + val;
		total_other += val;

		all_counts[i] = all_counts[i] + (val*3.0);
		total_all += (val*3.0);


		bg_alpha_counts[i] = bg_alpha_counts[i] + val;
		bg_total_alpha += val;

		bg_beta_counts[i] = bg_beta_counts[i] + val;
		bg_total_beta += val;

		bg_other_counts[i] = bg_other_counts[i] + val;
		bg_total_other += val;

		bg_all_counts[i] = bg_all_counts[i] + (val*3.0);
		bg_total_all += (val*3.0);
	}


	for (int i = 0 ; i < 100; i++){
		alpha_counts[i] = alpha_counts[i] / total_alpha;
		beta_counts[i] = beta_counts[i] / total_beta;
		other_counts[i] = other_counts[i] / total_other;
		all_counts[i] = all_counts[i] / total_all;

		bg_alpha_counts[i] = bg_alpha_counts[i] / bg_total_alpha;
		bg_beta_counts[i] = bg_beta_counts[i] / bg_total_beta;
		bg_other_counts[i] = bg_other_counts[i] / bg_total_other;
		bg_all_counts[i] = bg_all_counts[i] / bg_total_all;

		alpha_pot[i] = -std::log(alpha_counts[i]/bg_alpha_counts[i]);
		beta_pot[i] = -std::log(beta_counts[i]/bg_beta_counts[i]);
		other_pot[i] = -std::log(other_counts[i]/bg_other_counts[i]);
		all_pot[i]  = -std::log(all_counts[i]/bg_all_counts[i]);

	}
	return true;
}

bool ca_density::output_training_results(std::ostream& output){
	for (int i = 0 ; i < 100; i++){
		cout << "ca_density\t" << i << "\talpha\t" << alpha_pot[i] << "\t" << alpha_counts[i] << "\t" << (alpha_counts[i]*total_alpha) << "\n";
		cout << "ca_density\t" << i << "\tbeta\t"  << beta_pot[i]  << "\t" << beta_counts[i]  << "\t" << (beta_counts[i]*total_beta)   << "\n";
		cout << "ca_density\t" << i << "\tother\t" << other_pot[i] << "\t" << other_counts[i] << "\t" << (other_counts[i]*total_other) << "\n";
		cout << "ca_density\t" << i << "\tall\t"   << all_pot[i]   << "\t" << all_counts[i]   << "\t" << (all_counts[i]*total_all)     << "\n";
	}
	return true;
}

std::istream& ca_density::load_data( std::istream& input ){
	alpha_pot.clear();
	beta_pot.clear();
	other_pot.clear();
	all_pot.clear();

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
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 4 ){

				string paraName = SplitVec[0];
				trim(paraName);

				if ( paraName.compare("ca_density") == 0 ) {
					int count = lexical_cast<int>(SplitVec[1]);
					string type = lexical_cast<string>(SplitVec[2]);
					double val = lexical_cast<double>(SplitVec[3]);

					if (type.compare("alpha") == 0){
						alpha_pot[count] = val;
					}
					else if (type.compare("beta") == 0){
						beta_pot[count] = val;
					}
					else if (type.compare("other") == 0){
						other_pot[count] = val;
					}
					else if (type.compare("all") == 0){
						//do nothing
					}
					else {
	                    cout << "ca_density: ERROR - unknown type name: " << type << endl;
					}

				}
				else {
                    cout << "ca_density: ERROR - unknown parameter name: " << paraName << endl;
                    //fatalError = true;
                    //cout << paraValue << endl;
                }


			}
		}
	}

	return input;
}



}
}
}
}


