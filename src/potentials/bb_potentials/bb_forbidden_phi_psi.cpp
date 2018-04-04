/*
 * bb_forbidden_phi_psi.cpp
 *
 *  Created on: 1 Nov 2010
 *      Author: jmacdona
 */

#include "bb_forbidden_phi_psi.h"

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

potential_shared_ptr new_bb_forbidden_phi_psi(){
	potential_shared_ptr ptr(new bb_forbidden_phi_psi());
	return ptr;
}



bb_forbidden_phi_psi::bb_forbidden_phi_psi(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_forbidden_phi_psi"));

	overall_freq_map.resize(grid_divs*grid_divs,0);
	energies.resize(grid_divs*grid_divs,0);
	totalResidueCount = 0;

	phi_psi_increment = (2.0*PI) / static_cast<double>(grid_divs);

	forbidden_cutoff = 4.0;
	hot_spot_cutoff = 0;
}

const unsigned int bb_forbidden_phi_psi::grid_divs = 20;

bool bb_forbidden_phi_psi::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:bb_forbidden_phi_psi");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	//count_cutoff = PRODART::ENV::get_option_value<int>("potential:bb:bb_forbidden_phi_psi:count_cutoff");

	return true;
}

void bb_forbidden_phi_psi::getPhiPsiBin(const int phi_psi_sector, int& phi, int& psi) const{

	/*
	 * 	getPhiGridRef( griddivisions ) +
	 * (griddivisions * getPsiGridRef( griddivisions ));
	 */

	psi = phi_psi_sector / grid_divs;

	phi = (phi_psi_sector - (psi * grid_divs));

	return;


}

double bb_forbidden_phi_psi::get_lower_bound(int bin) const{
    return (static_cast<double>(bin) *  phi_psi_increment) - PI;
}

double bb_forbidden_phi_psi::get_grid_centre(int bin) const{
    return (static_cast<double>(bin) *  phi_psi_increment) - PI + (phi_psi_increment/2);
}



double bb_forbidden_phi_psi::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const_residue_shared_ptr res_0;

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();
	int phipsiSector = 0;


	for (int i = 0; i < resCount ; i++){

		res_0 = pose_->get_residue(i);


		phipsiSector = bb_meta_dat->get_phi_psi_sector(i,20);


		if (!res_0->is_terminal()
				) {

			total_energy += energies[phipsiSector];
		}
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_forbidden_phi_psi::get_energy_residue_loop(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map,
			const bool_vector& res_loop_mask) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const_residue_shared_ptr res_0;

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();
	int phipsiSector = 0;


	for (int i = 1; i < resCount-1 ; i++){
		if (res_loop_mask[i] == true
				|| res_loop_mask[i-1] == true
				|| res_loop_mask[i+1] == true){
			res_0 = pose_->get_residue(i);


			phipsiSector = bb_meta_dat->get_phi_psi_sector(i,20);


			if (!res_0->is_terminal()
			) {

				total_energy += energies[phipsiSector];
			}
		}
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_forbidden_phi_psi::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	return this->get_energy(pose_meta_, energies_map);
}

void bb_forbidden_phi_psi::get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		bool_vector& vec) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	const_residue_shared_ptr res_0;

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();
	int phipsiSector = 0;


	for (int i = 0; i < resCount ; i++){

		res_0 = pose_->get_residue(i);


		phipsiSector = bb_meta_dat->get_phi_psi_sector(i,20);


		if (!res_0->is_terminal()
				) {

			if (energies[phipsiSector] > hot_spot_cutoff){
				vec[i] = true;
			}
		}
	}


}


std::istream& bb_forbidden_phi_psi::load_data( std::istream& input ){

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
                     && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 9 ){
                        string paraName = SplitVec[0];
                        trim(paraName);
                        const int bin = lexical_cast<int>(SplitVec[8]);
                        const double energy = lexical_cast<double>(SplitVec[2]);

                        this->energies[bin] = energy;

                    }

            }

    }

	return input;

}




}
}
}
}



