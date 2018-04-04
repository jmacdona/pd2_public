/*
 * bb_pp_bias_correction.cpp
 *
 *  Created on: 1 Nov 2010
 *      Author: jmacdona
 */

#include "bb_pp_bias_correction.h"

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

bb_pp_bias_correction_shared_ptr new_bb_pp_bias_correction(){
	bb_pp_bias_correction_shared_ptr ptr(new bb_pp_bias_correction());
	return ptr;
}



bb_pp_bias_correction::bb_pp_bias_correction(): histo({grid_divs, grid_divs}, {-PI, -PI}, {PI, PI}, 1, 1, true, false){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_pp_bias_corr"));

	overall_freq_map.resize(grid_divs*grid_divs,0);
	energies.resize(grid_divs*grid_divs,0);
	totalResidueCount = 0;
	phi_psi_increment = (2.0*PI) / static_cast<double>(grid_divs);
	hot_spot_cutoff = 1.6;

}

const unsigned int bb_pp_bias_correction::grid_divs = 20;

bool bb_pp_bias_correction::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:bb_pp_bias_correction");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	//count_cutoff = PRODART::ENV::get_option_value<int>("potential:bb:bb_pp_bias_correction:count_cutoff");

	return true;
}

void bb_pp_bias_correction::getPhiPsiBin(const int phi_psi_sector, int& phi, int& psi) const{

	/*
	 * 	getPhiGridRef( griddivisions ) +
	 * (griddivisions * getPsiGridRef( griddivisions ));
	 */

	psi = phi_psi_sector / grid_divs;

	phi = (phi_psi_sector - (psi * grid_divs));

	return;


}

double bb_pp_bias_correction::get_lower_bound(int bin) const{
    return (static_cast<double>(bin) *  phi_psi_increment) - PI;
}

double bb_pp_bias_correction::get_grid_centre(int bin) const{
    return (static_cast<double>(bin) *  phi_psi_increment) - PI + (phi_psi_increment/2);
}



double bb_pp_bias_correction::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const_residue_shared_ptr res_0;

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();
	int phipsiSector = 0;


	for (int i = 0; i < resCount ; i++){

		res_0 = pose_->get_residue(i);


		



		if (!res_0->is_terminal()
				) {
			
			atom_shared_ptr4_tuple phi_atoms = bb_meta_dat->get_pose()->get_phi_atoms(i);
			atom_shared_ptr4_tuple psi_atoms = bb_meta_dat->get_pose()->get_psi_atoms(i);
			
			if (get<0>(phi_atoms)->isActiveAndSet() && get<1>(phi_atoms)->isActiveAndSet() &&
				get<2>(phi_atoms)->isActiveAndSet() && get<3>(phi_atoms)->isActiveAndSet() &&
				get<0>(psi_atoms)->isActiveAndSet() && get<1>(psi_atoms)->isActiveAndSet() &&
				get<2>(psi_atoms)->isActiveAndSet() && get<3>(psi_atoms)->isActiveAndSet()){
			
			//phipsiSector = bb_meta_dat->get_phi_psi_sector(i,20);
			double phi = bb_meta_dat->get_phi(i);
			double psi = bb_meta_dat->get_psi(i);
			
			double df_dx = 0, df_dy = 0;
			
			const double new_energy = histo.get_interp_val_no_recalc(phi,psi,df_dx, df_dy);
			//const double old_energy = energies[phipsiSector];
			//std::cout << "OLD_NEW:\t" << old_energy << "\t" << new_energy << std::endl;
			total_energy += new_energy;
			}


		}
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_pp_bias_correction::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();

	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const_residue_shared_ptr res_0;

	const int resCount = pose_->get_residue_count();//protein.getResidueCount();
	//int phipsiSector = 0;


	for (int i = 0; i < resCount ; i++){

		res_0 = pose_->get_residue(i);


		//phipsiSector = bb_meta_dat->get_phi_psi_sector(i,20);




		if (!res_0->is_terminal()
		) {

					atom_shared_ptr4_tuple phi_atoms = bb_meta_dat->get_pose()->get_phi_atoms(i);
					atom_shared_ptr4_tuple psi_atoms = bb_meta_dat->get_pose()->get_psi_atoms(i);
			
			if (get<0>(phi_atoms)->isActiveAndSet() && get<1>(phi_atoms)->isActiveAndSet() &&
				get<2>(phi_atoms)->isActiveAndSet() && get<3>(phi_atoms)->isActiveAndSet() &&
				get<0>(psi_atoms)->isActiveAndSet() && get<1>(psi_atoms)->isActiveAndSet() &&
				get<2>(psi_atoms)->isActiveAndSet() && get<3>(psi_atoms)->isActiveAndSet()){
			
					double phi = bb_meta_dat->get_phi(i);
					const vector3d phi_vec[4] = {get<0>(phi_atoms)->get_coords(), get<1>(phi_atoms)->get_coords(),
							get<2>(phi_atoms)->get_coords(), get<3>(phi_atoms)->get_coords()};
					const int phi_indices[4] = {get<0>(phi_atoms)->get_seq_num(), get<1>(phi_atoms)->get_seq_num(),
							get<2>(phi_atoms)->get_seq_num(), get<3>(phi_atoms)->get_seq_num()};
			
					double psi = bb_meta_dat->get_psi(i);
					const vector3d psi_vec[4] = {get<0>(psi_atoms)->get_coords(), get<1>(psi_atoms)->get_coords(),
									get<2>(psi_atoms)->get_coords(), get<3>(psi_atoms)->get_coords()};
					const int psi_indices[4] = {get<0>(psi_atoms)->get_seq_num(), get<1>(psi_atoms)->get_seq_num(),
									get<2>(psi_atoms)->get_seq_num(), get<3>(psi_atoms)->get_seq_num()};



					vector3d at_phi_grad_vecs[4];
					dihedral_p1_cosine_deriv(
							phi_vec[0],
							phi_vec[1],
							phi_vec[2],
							phi_vec[3],
							phi,
							at_phi_grad_vecs[0]);
					dihedral_p1_cosine_deriv(
							phi_vec[3],
							phi_vec[2],
							phi_vec[1],
							phi_vec[0],
							phi,
							at_phi_grad_vecs[3]);
					dihedral_p2_cosine_deriv(
							phi_vec[0],
							phi_vec[1],
							phi_vec[2],
							phi_vec[3],
							phi,
							at_phi_grad_vecs[1]);
					dihedral_p2_cosine_deriv(
							phi_vec[3],
							phi_vec[2],
							phi_vec[1],
							phi_vec[0],
							phi,
							at_phi_grad_vecs[2]);

					vector3d at_psi_grad_vecs[4];
					dihedral_p1_cosine_deriv(
							psi_vec[0],
							psi_vec[1],
							psi_vec[2],
							psi_vec[3],
							psi,
							at_psi_grad_vecs[0]);
					dihedral_p1_cosine_deriv(
							psi_vec[3],
							psi_vec[2],
							psi_vec[1],
							psi_vec[0],
							psi,
							at_psi_grad_vecs[3]);
					dihedral_p2_cosine_deriv(
							psi_vec[0],
							psi_vec[1],
							psi_vec[2],
							psi_vec[3],
							psi,
							at_psi_grad_vecs[1]);
					dihedral_p2_cosine_deriv(
							psi_vec[3],
							psi_vec[2],
							psi_vec[1],
							psi_vec[0],
							psi,
							at_psi_grad_vecs[2]);

				double df_dphi = 0, df_dpsi = 0;

			const double new_energy = histo.get_interp_val_no_recalc(phi,psi,df_dphi, df_dpsi);
			//const double old_energy = energies[phipsiSector];
			//std::cout << "OLD_NEW:\t" << old_energy << "\t" << new_energy << std::endl;
			total_energy += new_energy;

			for (int i = 0; i < 4; i++){
				grad[phi_indices[i]] += df_dphi * at_phi_grad_vecs[i];
				grad[psi_indices[i]] += df_dpsi * at_psi_grad_vecs[i];
				//cout << "ana: " << i << "\t" << dE_ddih * at_grad_vecs[i] << endl;
			}

			}
		}
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

void bb_pp_bias_correction::get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
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

			if ( energies[phipsiSector] > hot_spot_cutoff) {
				vec[i] = true;
			}
		}
	}


}


std::istream& bb_pp_bias_correction::load_data( std::istream& input ){


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
                     && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 3 ){
                        string paraName = SplitVec[0];
                        trim(paraName);
                        const int bin = lexical_cast<int>(SplitVec[1]);
                        const double energy = lexical_cast<double>(SplitVec[2]);
                        int phi_bin = 0, psi_bin = 0;
                        getPhiPsiBin(bin, phi_bin, psi_bin);

                        if ( paraName.compare("PPBC") == 0 ) {

                        	this->energies[bin] = energy;
                        	histo.set_count({(unsigned int)phi_bin, (unsigned int)psi_bin}, energy);

                        }



                    }

            }

    }

    histo.precal_bicubic_arr();

	return input;

}



bool bb_pp_bias_correction::addPseudoCounts( const double sumCount){
	double_vector::iterator iter;

	for (iter = this->overall_freq_map.begin(); iter != this->overall_freq_map.end(); iter++){
		*iter += sumCount;
		this->totalResidueCount += sumCount;
	}

	return true;
}


bool bb_pp_bias_correction::calculateScores( void ){

	const int numBins = overall_freq_map.size();
	energies.resize(numBins, 0);

	const double expected_count = totalResidueCount / static_cast<double>(numBins);

	for (int bin = 0; bin < numBins; bin++){
		this->energies[bin] = -std::log(overall_freq_map[bin] / expected_count);
	}

	return true;


}


void bb_pp_bias_correction::addToDB(const double phi, const double psi){
	if (phi >= - PI && phi <= PI
			&& psi >= - PI && psi <= PI){
		const int phipsiSector = pose_meta_interface::get_phi_psi_sector(phi, psi, grid_divs);
		if ((phipsiSector < (int)overall_freq_map.size())
				&& phipsiSector >=0 ){
			this->overall_freq_map[phipsiSector]++;
			totalResidueCount = totalResidueCount +1;
		}
		else {
			cerr << "ERROR: bin out of range" << endl;
			cerr << phipsiSector << "\t" << phi << "\t" << psi << "\t" << endl;
		}
	}
	else {
		cerr << "ERROR: dihedrals out of range" << endl;
		cerr << phi << "\t" << psi << "\t" << endl;
	}
}


std::ostream& bb_pp_bias_correction::output_phi_psi_info(std::ostream& output){
	/*
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	cout << "#phi_psi_bias_correction: " << asctime (timeinfo) << "#" << endl;
	 */



	const int numBins = energies.size();

	for (int bin = 0; bin < numBins; bin++){

		int phi = 0, psi = 0;

		getPhiPsiBin(bin, phi, psi );
		double lb_phi = get_lower_bound(phi);
		double lb_psi = get_lower_bound(psi);

		double c_phi = get_grid_centre(phi);
		double c_psi = get_grid_centre(psi);

		output << c_phi * (180.0 / PI) << "\t"
			 << c_psi * (180.0 / PI) << "\t"
			 << energies[bin] << "\t"
			 << overall_freq_map[bin] << "\t"
			 << lb_phi * (180.0 / PI) << "\t"
			 << lb_psi * (180.0 / PI) << "\t"
			 << phi << "\t"
			 << psi << "\t"
			 << bin << "\t"
		     << endl;

	}

	return output;
}



}
}
}
}

