/*
 * bb_frag3_mq.cpp
 *
 *  Created on: Oct 31, 2010
 *      Author: jmacdon
 */

#include "bb_frag3_mq.h"

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
namespace BB{

potential_shared_ptr new_bb_frag3_mq(){
	potential_shared_ptr ptr(new bb_frag3_mq());
	return ptr;
}



bb_frag3_mq::bb_frag3_mq(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_frag3_mq"));

	frag3_counts.clear();
	num_phi_psi_bins = 10;
	num_omega_bins = 5;
	count_cutoff = 2; // was 5

}

bool bb_frag3_mq::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:bb_frag3_mq");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	count_cutoff = PRODART::ENV::get_option_value<int>("potential:bb:bb_frag3_mq:count_cutoff");

	return true;
}

int bb_frag3_mq::get_dih_bin(const double dih, const int num_bins) const{
    const double increment = (2.0*PI) / static_cast<double>(num_bins);
    return static_cast<int>((dih + PI) / increment);
}




double bb_frag3_mq::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();

	for (int i = 1; i < resCount - 2; i++){
		const_residue_shared_ptr res_i = pose_->get_residue(i);
		const_residue_shared_ptr res_i_p1 = pose_->get_residue(i+1);
		const_residue_shared_ptr res_i_p2 = pose_->get_residue(i+2);

		if (!res_i->is_terminal()
				&& !res_i_p1->is_terminal()
				&& !res_i_p2->is_terminal()){

			const int phi_i = get_dih_bin(bb_meta_dat->get_phi(i), this->num_phi_psi_bins);
			const int psi_i = get_dih_bin(bb_meta_dat->get_psi(i), this->num_phi_psi_bins);

			const int omega_i_p1 = static_cast<int>(bb_meta_dat->is_trans(i+1));//get_dih_bin(bb_meta_dat->get_omega(i+1), this->num_omega_bins);
			const int phi_i_p1 = get_dih_bin(bb_meta_dat->get_phi(i+1), this->num_phi_psi_bins);
			const int psi_i_p1 = get_dih_bin(bb_meta_dat->get_psi(i+1), this->num_phi_psi_bins);

			const int omega_i_p2 = static_cast<int>(bb_meta_dat->is_trans(i+2));//get_dih_bin(bb_meta_dat->get_omega(i+2), this->num_omega_bins);
			const int phi_i_p2 = get_dih_bin(bb_meta_dat->get_phi(i+2), this->num_phi_psi_bins);
			const int psi_i_p2 = get_dih_bin(bb_meta_dat->get_psi(i+2), this->num_phi_psi_bins);


			stringstream keySStr(stringstream::in | stringstream::out);


			keySStr << phi_i << "-"
					<< psi_i << "-"
					<< omega_i_p1 << "-"
					<< phi_i_p1 << "-"
					<< psi_i_p1 << "-"
					<< omega_i_p2 << "-"
					<< phi_i_p2 << "-"
					<< psi_i_p2;// << "-";


			if ((frag3_counts.find(keySStr.str()) != frag3_counts.end() ? frag3_counts.find(keySStr.str())->second : 0) < count_cutoff ){		//frag3_counts[keySStr.str()] < count_cutoff){
				total_energy++;
			}



		}

	}


	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_frag3_mq::get_energy_residue_loop(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map,
			const bool_vector& res_loop_mask) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();

	for (int i = 1; i < resCount - 2; i++){
		if (res_loop_mask[i+1]){
			const_residue_shared_ptr res_i = pose_->get_residue(i);
			const_residue_shared_ptr res_i_p1 = pose_->get_residue(i+1);
			const_residue_shared_ptr res_i_p2 = pose_->get_residue(i+2);

			if (!res_i->is_terminal()
					&& !res_i_p1->is_terminal()
					&& !res_i_p2->is_terminal()){

				const int phi_i = get_dih_bin(bb_meta_dat->get_phi(i), this->num_phi_psi_bins);
				const int psi_i = get_dih_bin(bb_meta_dat->get_psi(i), this->num_phi_psi_bins);

				const int omega_i_p1 = static_cast<int>(bb_meta_dat->is_trans(i+1));//get_dih_bin(bb_meta_dat->get_omega(i+1), this->num_omega_bins);
				const int phi_i_p1 = get_dih_bin(bb_meta_dat->get_phi(i+1), this->num_phi_psi_bins);
				const int psi_i_p1 = get_dih_bin(bb_meta_dat->get_psi(i+1), this->num_phi_psi_bins);

				const int omega_i_p2 = static_cast<int>(bb_meta_dat->is_trans(i+2));//get_dih_bin(bb_meta_dat->get_omega(i+2), this->num_omega_bins);
				const int phi_i_p2 = get_dih_bin(bb_meta_dat->get_phi(i+2), this->num_phi_psi_bins);
				const int psi_i_p2 = get_dih_bin(bb_meta_dat->get_psi(i+2), this->num_phi_psi_bins);


				stringstream keySStr(stringstream::in | stringstream::out);


				keySStr << phi_i << "-"
						<< psi_i << "-"
						<< omega_i_p1 << "-"
						<< phi_i_p1 << "-"
						<< psi_i_p1 << "-"
						<< omega_i_p2 << "-"
						<< phi_i_p2 << "-"
						<< psi_i_p2;// << "-";


				if ((frag3_counts.find(keySStr.str()) != frag3_counts.end() ? frag3_counts.find(keySStr.str())->second : 0) < count_cutoff ){		//frag3_counts[keySStr.str()] < count_cutoff){
					total_energy++;
				}



			}
		}

	}


	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_frag3_mq::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	return this->get_energy(pose_meta_, energies_map);
}

void bb_frag3_mq::get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		bool_vector& vec) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();

	for (int i = 1; i < resCount - 2; i++){
		const_residue_shared_ptr res_i = pose_->get_residue(i);
		const_residue_shared_ptr res_i_p1 = pose_->get_residue(i+1);
		const_residue_shared_ptr res_i_p2 = pose_->get_residue(i+2);

		if (!res_i->is_terminal()
				&& !res_i_p1->is_terminal()
				&& !res_i_p2->is_terminal()){

			const int phi_i = get_dih_bin(bb_meta_dat->get_phi(i), this->num_phi_psi_bins);
			const int psi_i = get_dih_bin(bb_meta_dat->get_psi(i), this->num_phi_psi_bins);

			const int omega_i_p1 = static_cast<int>(bb_meta_dat->is_trans(i+1));//get_dih_bin(bb_meta_dat->get_omega(i+1), this->num_omega_bins);
			const int phi_i_p1 = get_dih_bin(bb_meta_dat->get_phi(i+1), this->num_phi_psi_bins);
			const int psi_i_p1 = get_dih_bin(bb_meta_dat->get_psi(i+1), this->num_phi_psi_bins);

			const int omega_i_p2 = static_cast<int>(bb_meta_dat->is_trans(i+2));//get_dih_bin(bb_meta_dat->get_omega(i+2), this->num_omega_bins);
			const int phi_i_p2 = get_dih_bin(bb_meta_dat->get_phi(i+2), this->num_phi_psi_bins);
			const int psi_i_p2 = get_dih_bin(bb_meta_dat->get_psi(i+2), this->num_phi_psi_bins);


			stringstream keySStr(stringstream::in | stringstream::out);


			keySStr << phi_i << "-"
					<< psi_i << "-"
					<< omega_i_p1 << "-"
					<< phi_i_p1 << "-"
					<< psi_i_p1 << "-"
					<< omega_i_p2 << "-"
					<< phi_i_p2 << "-"
					<< psi_i_p2;// << "-";


			if ((frag3_counts.find(keySStr.str()) != frag3_counts.end() ? frag3_counts.find(keySStr.str())->second : 0) < count_cutoff ){		//frag3_counts[keySStr.str()] < count_cutoff){
				//vec[i] = true;
				vec[i+1] = true;
				//vec[i+2] = true;
			}



		}

	}


}


std::istream& bb_frag3_mq::load_data( std::istream& input ){

	string lineStr;


	long length, lineNum = 0 ;


	string_vector SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 2 ){

				string key = SplitVec[0];
				trim(key);

				const double count = lexical_cast<double>(SplitVec[1]);

				this->frag3_counts[key] = count;

			}
		}
	}

	return input;

}



}
}
}
}

