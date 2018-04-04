/*
 * misc_protocols.cpp
 *
 *  Created on: 18 Nov 2010
 *      Author: jmacdona
 */

#include "misc_protocols.h"
#include "protocols.h"

using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::POTENTIALS::BB;
using namespace PRODART::POSE::SIM;
using namespace PRODART::ROTAMERS;

using boost::split;
using boost::filesystem::path;
using std::set;
//using namespace boost::filesystem;

namespace PRODART {
namespace PROTOCOLS{
namespace MISC{

void output_rosetta_contraints_file_polyALA(PRODART::POSE::pose_shared_ptr protein, std::ostream& output){

	protein->index();

	bb_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_bb_pose_meta(protein);
	PRODART::POSE::three_state_sec_struct_vector secs = pose_meta_->get_sec_struct();

	const int resCount = protein->get_residue_count();

	for (int i = 0 ; i< resCount; i++){
		for (int j = i + 5 ; j < resCount; j++){

			if ( (secs[i] == ss3_HELIX && secs[j] == ss3_STRAND)
					|| (secs[j] == ss3_HELIX && secs[i] == ss3_STRAND) ){
				const double dist = (protein->get_bb_atom(POSE::CA, i)->get_coords() - protein->get_bb_atom(POSE::CA, j)->get_coords()).mod();

				if (dist < ENV::get_option_value<double>("rosetta_contraints:cst_ca_cutoff")){
					int res_i = protein->get_residue(i)->get_internal_residue_index() + 1;//get_pdb_residue_index();
					int res_j = protein->get_residue(j)->get_internal_residue_index() + 1;//get_pdb_residue_index();

					//trim(res_i);
					//trim(res_j);

					output << "AtomPair CA " << res_i << " CA " << res_j << " HARMONIC " << dist << " " << "1" << "\n";
				}

			}


		}
	}

}

void output_rosetta_contraints_file_all(PRODART::POSE::pose_shared_ptr protein, std::ostream& output){

	protein->index();

	bb_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_bb_pose_meta(protein);
	//PRODART::POSE::three_state_sec_struct_vector secs = pose_meta_->get_sec_struct();

	const int resCount = protein->get_residue_count();

	for (int i = 0 ; i< resCount; i++){
		for (int j = i + 5 ; j < resCount; j++){


			const double dist = (protein->get_bb_atom(POSE::CA, i)->get_coords() - protein->get_bb_atom(POSE::CA, j)->get_coords()).mod();

			if (dist < ENV::get_option_value<double>("rosetta_contraints:cst_ca_cutoff")){
				int res_i = protein->get_residue(i)->get_internal_residue_index() + 1;//get_pdb_residue_index();
				int res_j = protein->get_residue(j)->get_internal_residue_index() + 1;//get_pdb_residue_index();

				//trim(res_i);
				//trim(res_j);

				output << "AtomPair CA " << res_i << " CA " << res_j << " HARMONIC " << dist << " " << "1" << "\n";
			}




		}
	}
}

void output_rosetta_contraints_file_all(PRODART::POSE::pose_shared_ptr protein, bool_vector loop_mask, std::ostream& output){
	protein->index();

	bb_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_bb_pose_meta(protein);
	//PRODART::POSE::three_state_sec_struct_vector secs = pose_meta_->get_sec_struct();

	const int resCount = protein->get_residue_count();

	for (int i = 0 ; i< resCount; i++){
		for (int j = i + 5 ; j < resCount; j++){


			if (loop_mask[i] && loop_mask[j]){
				const double dist = (protein->get_bb_atom(POSE::CA, i)->get_coords() - protein->get_bb_atom(POSE::CA, j)->get_coords()).mod();

				if (dist < ENV::get_option_value<double>("rosetta_contraints:cst_ca_cutoff")){
					int res_i = protein->get_residue(i)->get_internal_residue_index() + 1;//get_pdb_residue_index();
					int res_j = protein->get_residue(j)->get_internal_residue_index() + 1;//get_pdb_residue_index();

					//trim(res_i);
					//trim(res_j);

					output << "AtomPair CA " << res_i << " CA " << res_j << " HARMONIC " << dist << " " << "1" << "\n";
				}
			}




		}
	}
}

void output_rosetta_sc_sc_contraints_file(PRODART::POSE::pose_shared_ptr protein, bool_vector loop_mask, std::ostream& output){
	protein->index();

	bb_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_bb_pose_meta(protein);
	//PRODART::POSE::three_state_sec_struct_vector secs = pose_meta_->get_sec_struct();

	const int resCount = protein->get_residue_count();

	for (int i = 0 ; i< resCount; i++){
		for (int j = i + 5 ; j < resCount; j++){


			if (loop_mask[i] && loop_mask[j]){

				const int res_i = protein->get_residue(i)->get_internal_residue_index() + 1;//get_pdb_residue_index();
				const int res_j = protein->get_residue(j)->get_internal_residue_index() + 1;//get_pdb_residue_index();

				sidechain_shared_ptr sc_i = protein->get_residue(i)->get_sidechain();
				sidechain_shared_ptr sc_j = protein->get_residue(j)->get_sidechain();

				const int sc_i_at_count = sc_i->get_atom_count();
				const int sc_j_at_count = sc_j->get_atom_count();

				for (int sc_i_num = 0; sc_i_num < sc_i_at_count; sc_i_num++){
					for (int sc_j_num = 0; sc_j_num < sc_j_at_count; sc_j_num++){
						atom_shared_ptr at_i = sc_i->get_atom(sc_i_num);
						atom_shared_ptr at_j = sc_j->get_atom(sc_j_num);
						if (at_i->isActiveAndSet() && at_j->isActiveAndSet()){
							const double dist = (at_i->get_coords() - at_j->get_coords()).mod();
							if (dist < ENV::get_option_value<double>("rosetta_contraints:cst_ca_cutoff")){
								output << "AtomPair " << at_i->get_type().get_label() << " " << res_i
										<< " " << at_j->get_type().get_label() << res_j << " HARMONIC " << dist << " " << "10" << "\n";
								//output << "AtomPair CA " << res_i << " CA " << res_j << " HARMONIC " << dist << " " << "2" << "\n";
							}
						}
					}
				}

			}




		}
	}
}

void trim_terminal_tails(PRODART::POSE::pose_shared_ptr protein, const int max_tail){

	protein->index();

	const int resCount = protein->get_residue_count();


	for (int i = 0; i < resCount; i++){
		atom_shared_ptr at = protein->get_bb_atom(POSE::H, i);
		at->set_type(atom_type("H"));
	}

	bb_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_bb_pose_meta(protein);
	PRODART::POSE::three_state_sec_struct_vector secs = pose_meta_->get_sec_struct();


	int c_term_del = 0, n_term_del =0;



	for (int i = resCount-1; i != 0; i--){
		if (secs[i] == PRODART::POSE::ss3_OTHER ){
			bool ok = true;
			for (int j = i; j >= i - max_tail; j--){
				if (secs[j] != ss3_OTHER){
					ok = false;
				}
			}
			if (ok){
				c_term_del++;
			}
		}
		else {
			break;
		}
	}

	for (int i = 0; i < resCount; i++){
		if (secs[i] == PRODART::POSE::ss3_OTHER ){
			bool ok = true;
			for (int j = i; j <= i + max_tail; j++){
				if (secs[j] != ss3_OTHER){
					ok = false;
				}
			}
			if (ok){
				n_term_del++;
			}
		}
		else {
			break;
		}
	}

	for (int i = resCount-1; i>resCount-1-c_term_del; i--){
		protein->delete_residue(i);
	}

	for (int i = 0; i<n_term_del; i++){
		protein->delete_residue(0);
	}

	cout << "deleted " << n_term_del << " n-term residues and";
	cout << " " << c_term_del << " c-term residues" << endl;


	protein->index();

}

void inactivate_unset_GLY_CBs(PRODART::POSE::pose_shared_ptr protein){

	const int rescount = protein->get_residue_count();

	for (int i = 0; i < rescount; i++){
		if (protein->get_residue(i)->get_type() == residue_type("GLY")){
			protein->get_bb_atom(CB, i)->setActive(false);
			protein->get_bb_atom(CB, i)->setSet(false);
		}
	}


}

void clean_pose(PRODART::POSE::pose_shared_ptr protein, bool remove_bb_H){
	const int at_count = protein->get_all_atom_count();
	for (int i = 0; i < at_count; i++ ){
		atom_shared_ptr atm = protein->get_atom(i);
		if (atm){
			atm->set_occupancy(1.00);
			atm->set_b_factor(20.00);
			if (remove_bb_H && atm->get_type() == atom_type("H")){
				atm->setActive(false);
				atm->setSet(false);
			}
		}
	}
	protein->index();

	const int res_count = protein->get_residue_count();
	for (int i = 0; i < res_count; i++ ){
		residue_shared_ptr res = protein->get_residue(i);
		if (res){
			if (!residue_type_map::Instance()->is_standard_type(res->get_type())){
				res->set_type(residue_type("ALA"));
			}
		}
	}

}

void to_bb_only(PRODART::POSE::pose_shared_ptr protein){

	const int rescount = protein->get_residue_count();
	for (int i = 0 ; i < rescount; i++){
		protein->get_residue(i)->get_sidechain()->clear();
	}

}


void output_rama_histogram_from_raw_phi_psi( std::istream& input, std::ostream& output,
		const bool trans,
		const bool pro){
    string lineStr;

    bb_pp_bias_correction_shared_ptr db = new_bb_pp_bias_correction();
	cout << "adding pseudocounts" << endl;
    db->addPseudoCounts(0.1);

    unsigned long length, lineNum = 0 ;


    string_vector SplitVec;
	cout << "reading file..." << endl;
    while ( !input.eof() ) {
            getline(input, lineStr);
            string resStr;
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {
                    split( SplitVec, lineStr, is_any_of("\t") );
                    if ( SplitVec[0].substr(0,1).compare("#") != 0
                     && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 5 ){

                    	string resname = SplitVec[1];
                    	trim(resname);
                        const double phi = degrees_to_radians(lexical_cast<double>(SplitVec[2]));
                        const double psi = degrees_to_radians(lexical_cast<double>(SplitVec[3]));
                        const double omega = degrees_to_radians(lexical_cast<double>(SplitVec[4]));

                        POSE::four_state_sec_struct ss =  get_phi_psi_omega_sector( phi, psi,  omega);




                        if (ss != ss4_CIS
                        		&& trans
                        		&& ((resname.compare("PRO") != 0) || pro )){
                            //cout << lineNum << "\t" << flush;
                        	db->addToDB(phi, psi);
                        	//PRINT_EXPR("a");
                        }
                        else if (ss == ss4_CIS
                        		&& !trans
                        		&& (resname.compare("PRO") == 0) ){
                        	db->addToDB(phi, psi);
                        	//PRINT_EXPR("b");
                        }
                        else {
                        	//PRINT_EXPR(ss);
                        	//PRINT_EXPR(trans);
                        	//PRINT_EXPR(resname);
                        	//PRINT_EXPR("c");
                        }



                    }

            }

    }

	cout << "calculating scores..." << endl;
    db->calculateScores();

	cout << "outputting histogram..." << endl;
    db->output_phi_psi_info(output);



}

void output_rama_histogram_from_raw_phi_psi_filelist( std::istream& list_input, std::ostream& output,
		const bool trans){
    string lineStr;

    bb_pp_bias_correction_shared_ptr db = new_bb_pp_bias_correction();
	cout << "adding pseudocounts" << endl;
    db->addPseudoCounts(0.1);

    unsigned long list_length, list_lineNum = 0 ;


    string_vector list_SplitVec;
	cout << "reading filelist..." << endl;
    while ( !list_input.eof() ) {
            getline(list_input, lineStr);
            string resStr;
            list_lineNum++;

            list_length = lineStr.length();

            //cout << endl << list_lineNum << " " << list_length << " ";

            if (list_length > 0) {
                    split( list_SplitVec, lineStr, is_any_of("\t") );
                    if ( list_SplitVec[0].substr(0,1).compare("#") != 0
                     && list_SplitVec[0].substr(0,1).compare("") != 0 && list_SplitVec.size() >= 1 ){

                    	string filename = list_SplitVec[0];
                    	//cout << "opening file " << filename << endl;
        				ifstream input(filename.c_str(), ios::in);

        				if (input.is_open()){
        				    string lineStr;
        				    unsigned long length, lineNum = 0 ;
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
        				                     && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 5 ){

        				                    	string resname = SplitVec[1];
        				                    	trim(resname);

        				                        const double phi = degrees_to_radians(lexical_cast<double>(SplitVec[2]));
        				                        const double psi = degrees_to_radians(lexical_cast<double>(SplitVec[3]));
        				                        const double omega = degrees_to_radians(lexical_cast<double>(SplitVec[4]));

        				                        POSE::four_state_sec_struct ss =  get_phi_psi_omega_sector( phi, psi,  omega);




        				                        if (ss != ss4_CIS
        				                        		&& trans
        				                        		&& resname.compare("PRO") != 0){
        				                            //cout << lineNum << "\t" << flush;
        				                        	db->addToDB(phi, psi);
        				                        }
        				                        else if (ss == ss4_CIS
        				                        		&& !trans
        				                        		&& resname.compare("PRO") == 0){
        				                        	db->addToDB(phi, psi);
        				                        }



        				                    }

        				            }

        				    }

        				}
        				else {
        					cerr << "output_rama_histogram_from_raw_phi_psi_filelist: ERROR: can't open file: " << filename << endl;
        				}

                    }

            }

    }

	cout << "calculating scores..." << endl;
    db->calculateScores();

	cout << "outputting histogram..." << endl;
    db->output_phi_psi_info(output);

}


void output_raw_phi_psi_omega(PRODART::POSE::pose_shared_ptr protein, std::ostream& output, std::string label){

	string lineStr;
	unsigned long length, lineNum = 0 ;
	string_vector SplitVec;




	//PRODART::POSE::pose_shared_ptr protein = PRODART::POSE::new_pose();
	//protein->loadPdb(in_pdb);
	bb_pose_meta_shared_ptr meta_ = new_bb_pose_meta(protein);

	const int rescount = protein->get_residue_count();

	for (int i = 0 ; i < rescount; i++){

		residue_shared_ptr res  = protein->get_residue(i);
		chain_shared_ptr ch = res->get_chain();

		if (!res->is_terminal() && ch->isPeptide()){

			const double phi = meta_->get_phi(i);
			const double psi = meta_->get_psi(i);
			const double omega = meta_->get_omega(i);

			output << res->get_trimmed_pdb_residue_index() << "\t"
					<< res->get_chain()->getChainID() << "\t"
					<< res->get_type().get_label3() << "\t"
					<< radians_to_degrees(phi) << "\t"
					<< radians_to_degrees(psi) << "\t"
					<< radians_to_degrees(omega) << "\t"
					<< "middle\t"
					<< label
					<< endl;

		}

	}



}


void output_raw_phi_psi_omega_from_list(std::istream& input, std::ostream& output){

    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;

	cout << "reading file..." << endl;
    while ( !input.eof() ) {
            getline(input, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

            	ifstream in_pdb(lineStr.c_str(), ios::in);

            	if (in_pdb.is_open()){

            		PRODART::POSE::pose_shared_ptr protein = PRODART::POSE::new_pose();
            		protein->loadPdb(in_pdb);

            		output_raw_phi_psi_omega(protein, output, lineStr);

            	} //


            }

    }

}

std::istream& next_ENDMDL(std::istream& input){

	string lineStr;
	string record;


	//int lineLength = 0;

	while ( !input.eof() ) {
		getline(input, lineStr);

		//lineLength = lineStr.length();
		record = lineStr.substr(0,6);

		if (record.compare( "ENDMDL" ) == 0 ) {

			return input;
		}
	}
	return input;
}

PRODART::POSE::pose_shared_ptr get_single_mdl(const int mdl_num, std::istream& input){

	pose_shared_ptr rtn_pose = new_pose();

	bool got_it = false;

	int count = 1;

	while ( !input.eof() && !got_it){



		if (count == mdl_num){
			pose_shared_ptr pose_ = new_pose();
			pose_->loadPdb(input);
			rtn_pose = pose_;
			got_it = true;
		}
		else {
			next_ENDMDL(input);
		}

		count++;
	}

	return rtn_pose;
}


class verified_load_pose_exception : public std::exception{
  virtual const char* what() const throw()
  {
    return "ERROR: PDB file contained errors";
  }
};


PRODART::POSE::pose_shared_ptr unverified_load_pose( const std::string filename ){
	PRODART::POSE::pose_shared_ptr test_pose = PRODART::POSE::new_pose();
	std::ifstream protein_file(filename.c_str(), ios::in);

	bool loadOK = true;

	if (protein_file.is_open()){
		test_pose->loadPdb(protein_file);
	}
	else {
		cerr << "unverified_load_pose: ERROR: can't open PDB file: " << filename
		                                           << endl;
		loadOK = false;
	}
	protein_file.close();

	if (loadOK == false){
		throw verified_load_pose_exception();
	}

	return test_pose;
}

PRODART::POSE::pose_shared_ptr verified_load_ca_pose( const std::string filename ){
	PRODART::POSE::pose_shared_ptr test_pose = PRODART::POSE::new_pose();
	std::ifstream protein_file(filename.c_str(), ios::in);

	bool loadOK = true;

	if (protein_file.is_open()){
		test_pose->loadPdb(protein_file);
		loadOK = validate_ca_pose(test_pose);
		if (!loadOK){
			cerr << "verified_load_ca_pose: ERROR: PDB file can not be verified: " << filename << endl;
		}
	}
	else {
		cerr << "verified_load_ca_pose: ERROR: can't open PDB file: " << filename
		                                           << endl;
		loadOK = false;
	}
	protein_file.close();

	if (loadOK == false){
		throw verified_load_pose_exception();
	}

	return test_pose;
}

PRODART::POSE::pose_shared_ptr verified_load_bb_pose( const std::string filename ){
	PRODART::POSE::pose_shared_ptr test_pose = PRODART::POSE::new_pose();
	std::ifstream protein_file(filename.c_str(), ios::in);

	bool loadOK = true;

	if (protein_file.is_open()){
		test_pose->loadPdb(protein_file);
		loadOK = validate_bb_pose(test_pose);
		if (!loadOK){
			cerr << "verified_load_bb_pose: ERROR: PDB file can not be verified: " << filename << endl;
		}
	}
	else {
		cerr << "verified_load_bb_pose: ERROR: can't open PDB file: " << filename
		                                           << endl;
		loadOK = false;
	}
	protein_file.close();

	if (loadOK == false){
		throw verified_load_pose_exception();
	}

	return test_pose;
}

class read_3state_secs_exception : public std::exception{
  virtual const char* what() const throw()
  {
    return "read_3state_secs_exception";
  }
};

PRODART::POSE::three_state_sec_struct_vector read_3state_secs(const std::string filename ){
	std::ifstream input(filename.c_str(), ios::in);
	PRODART::POSE::three_state_sec_struct_vector rtn_vec;
	rtn_vec.clear();

	bool loadOK = true;

	if (input.is_open()){
		input >> rtn_vec;
	}
	else {
		cerr << "read_3state_secs: ERROR: can't open sec struct file: " << filename
		                                           << endl;
		loadOK = false;
	}
	input.close();

	if (loadOK == false){
		throw read_3state_secs_exception();
	}

	return rtn_vec;
}

bool has_rosetta_match_atom_type(string res3, string atom_name ){
	if ( (res3.compare("ASP") == 0 && (atom_name.compare("OD2") == 0 || atom_name.compare("OD1") == 0  ))
			|| (res3.compare("GLU") == 0 && (atom_name.compare("OE2") == 0 || atom_name.compare("OE1") == 0  ))){
		// OOC DE
		return true;
	}
	else if (res3.compare("HIS") == 0 && (atom_name.compare("NE2") == 0 || atom_name.compare("ND1") == 0  )){
		// Nhis H
		return true;
	}
	else if (res3.compare("ARG") == 0 && (atom_name.compare("NH1") == 0 || atom_name.compare("NH2") == 0  )){
		// Narg R
		return true;
	}
	else if (res3.compare("CYS") == 0 && (atom_name.compare("SG") == 0 )){
		// S C
		return true;
	}
	return false;
}

string get_rosetta_match_atom_type(string res3, string atom_name ){
	if ( (res3.compare("ASP") == 0 && (atom_name.compare("OD2") == 0 || atom_name.compare("OD1") == 0  ))
			|| (res3.compare("GLU") == 0 && (atom_name.compare("OE2") == 0 || atom_name.compare("OE1") == 0  ))){
		// OOC DE
		return string("OOC");
	}
	else if (res3.compare("HIS") == 0 && (atom_name.compare("NE2") == 0 || atom_name.compare("ND1") == 0  )){
		// Nhis H
		return string("Nhis");
	}
	else if (res3.compare("ARG") == 0 && (atom_name.compare("NH1") == 0 || atom_name.compare("NH2") == 0  )){
		// Narg R
		return string("Narg");
	}
	else if (res3.compare("CYS") == 0 && (atom_name.compare("SG") == 0 )){
		// S C
		return string("S");
	}
	return string("ERROR!!!!!!!!!!!!!!!!!!!!");
}

string get_rosetta_match_residue1(string res3, string atom_name ){
	if ( (res3.compare("ASP") == 0 && (atom_name.compare("OD2") == 0 || atom_name.compare("OD1") == 0  ))
			|| (res3.compare("GLU") == 0 && (atom_name.compare("OE2") == 0 || atom_name.compare("OE1") == 0  ))){
		// OOC DE
		return string("DE");
	}
	else if (res3.compare("HIS") == 0 && (atom_name.compare("NE2") == 0 || atom_name.compare("ND1") == 0  )){
		// Nhis H
		return string("H");
	}
	else if (res3.compare("ARG") == 0 && (atom_name.compare("NH1") == 0 || atom_name.compare("NH2") == 0  )){
		// Narg R
		return string("R");
	}
	else if (res3.compare("CYS") == 0 && (atom_name.compare("SG") == 0 )){
		// S C
		return string("C");
	}
	return string("!!!ERROR!!!!!!!!!!!!!!!!!!!!");
}


void output_rosetta_match_csts(PRODART::POSE::atom_shared_ptr res1atm3,
		PRODART::POSE::atom_shared_ptr res1atm2,
		PRODART::POSE::atom_shared_ptr res1atm1,
		PRODART::POSE::atom_shared_ptr res2atm1,
		PRODART::POSE::atom_shared_ptr res2atm2,
		PRODART::POSE::atom_shared_ptr res2atm3,
		bool use_atom_types){

	vector3d r1a3 = res1atm3->get_coords();
	vector3d r1a2 = res1atm2->get_coords();
	vector3d r1a1 = res1atm1->get_coords();

	vector3d r2a1 = res2atm1->get_coords();
	vector3d r2a2 = res2atm2->get_coords();
	vector3d r2a3 = res2atm3->get_coords();

	const double distanceAB = (r1a1 - r2a1).mod();
	const double angle_A = angle(r1a2, r1a1, r2a1);
	const double angle_B = angle(r1a1, r2a1, r2a2);
	const double torsion_A = dihedral(r1a3, r1a2, r1a1, r2a1);
	const double torsion_AB = dihedral(r1a2, r1a1, r2a1, r2a2);
	const double torsion_B = dihedral(r1a1, r2a1, r2a2, r2a3);

	residue_shared_ptr res2 = res2atm1->get_residue();

	cout << "# block for " << res2->get_type().get_label3() << " "
			<< res2->get_trimmed_pdb_residue_index() << " "
			<< res2->get_chain()->getChainID() << "\n";
	cout << "# internal res number (for rosetta): " << res2->get_internal_residue_index() + 1 << "\n";
	cout << endl;

	cout << "CST::BEGIN\n";
	cout <<  "  TEMPLATE::   ATOM_MAP: 1 atom_name: " << res1atm1->get_type().get_trimmed_label() << " "
			<< res1atm2->get_type().get_trimmed_label() << " "
			<< res1atm3->get_type().get_trimmed_label() << " "
			<< endl;
	cout << "  TEMPLATE::   ATOM_MAP: 1 residue3: " << res1atm1->get_residue()->get_type().get_label3() << endl;

	cout << endl;

	if (use_atom_types == false || has_rosetta_match_atom_type(res2atm1->get_residue()->get_type().get_label3(), res2atm1->get_type().get_trimmed_label()) == false){
		cout <<  "  TEMPLATE::   ATOM_MAP: 2 atom_name: " << res2atm1->get_type().get_trimmed_label() << " "
				<< res2atm2->get_type().get_trimmed_label() << " "
				<< res2atm3->get_type().get_trimmed_label() << " "
				<< endl;
		cout << "  TEMPLATE::   ATOM_MAP: 2 residue3: " << res2atm1->get_residue()->get_type().get_label3() << endl;
	}
	else {
		cout <<  "  TEMPLATE::   ATOM_MAP: 2 atom_type: " << get_rosetta_match_atom_type(res2atm1->get_residue()->get_type().get_label3(), res2atm1->get_type().get_trimmed_label())
				<< endl;
		cout << "  TEMPLATE::   ATOM_MAP: 2 residue1: " << get_rosetta_match_residue1(res2atm1->get_residue()->get_type().get_label3(), res2atm1->get_type().get_trimmed_label()) << endl;
	}

	cout << endl;

	 cout <<  "  CONSTRAINT:: distanceAB:  " << distanceAB << " " << "0.1" << " " << "100.0" << " " << "1" << " " << "1" << endl;
	 cout <<  "  CONSTRAINT::    angle_A:  " << radians_to_degrees(angle_A) << " " << "5.0" << " " << "60.0" << " " << "360.0" << " " << "1" << endl;
	 cout <<  "  CONSTRAINT::    angle_B:  " << radians_to_degrees(angle_B) << " " << "5.0" << " " << "60.0" << " " << "360.0" << " " << "1" << endl;
	 cout <<  "  CONSTRAINT::  torsion_A:  " << radians_to_degrees(torsion_A) << " " << "10.0" << " " << "60.0" << " " << "360.0" << " " << "1" << endl;
	 cout <<  "  CONSTRAINT::  torsion_AB:  " << radians_to_degrees(torsion_AB) << " " << "10.0" << " " << "60.0" << " " << "360.0" << " " << "1" << endl;
	 cout <<  "  CONSTRAINT::  torsion_B:  " << radians_to_degrees(torsion_B) << " " << "10.0" << " " << "60.0" << " " << "360.0" << " " << "1" << endl;

	 cout << "CST::END\n";
	 cout << endl;
}

void output_rosetta_match_csts(PRODART::POSE::pose_shared_ptr protein,
		const std::string filename,
		bool use_atom_types ){

	ifstream input(filename.c_str(), ios::in);

	if (!input.is_open()){
		cerr << "ERROR: output_rosetta_match_csts: failed to open file: " << filename << endl;
		return;
	}

	PRODART::POSE::atom_shared_ptr res1atm3;
	PRODART::POSE::atom_shared_ptr res1atm2;
	PRODART::POSE::atom_shared_ptr res1atm1;
	PRODART::POSE::atom_shared_ptr res2atm1;
	PRODART::POSE::atom_shared_ptr res2atm2;
	PRODART::POSE::atom_shared_ptr res2atm3;

	string lineStr;

	unsigned long length, lineNum = 0 ;


	string_vector SplitVec;
	//cout << "reading file..." << endl;
	while ( !input.eof() ) {
		getline(input, lineStr);

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of(" \t") );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
					&& SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 3 ){

				string resid = SplitVec[0];
				char chainid = lexical_cast<char>(SplitVec[1]);
				atom_type atomid = atom_type(SplitVec[2]);
				atom_shared_ptr atm = protein->get_atom(atomid, resid, chainid);
				if (!atm){
					cerr << "ERROR: output_rosetta_match_csts: could not find atom: " << lineStr << endl;
					return;
				}
				lineNum++;

				switch (lineNum){
				case 1:
					res1atm3 = atm;
					break;
				case 2:
					res1atm2 = atm;
					break;
				case 3:
					res1atm1 = atm;
					break;
				case 4:
					res2atm1 = atm;
					break;
				case 5:
					res2atm2 = atm;
					break;
				case 6:
					res2atm3 = atm;
					break;
				default:
					cerr << "ERROR: output_rosetta_match_csts: only need 6 atoms defined: " << lineStr << endl;
					return;
				}


			}

		}

	}

	if (res1atm3 && res1atm2 && res1atm1 && res2atm1 && res2atm2 && res2atm3){
		output_rosetta_match_csts(res1atm3, res1atm2, res1atm1, res2atm1, res2atm2, res2atm3, use_atom_types);
	}
	else {
		cerr << "ERROR: output_rosetta_match_csts: need all 6 atoms defined, some are missing" << endl;
	}

}


void graft_motif_residues_to_scaffold(PRODART::POSE::pose_shared_ptr scaffold, PRODART::POSE::pose_shared_ptr in_motif){
	PRODART::POSE::pose_shared_ptr motif = in_motif->clone();
	const int mot_rescount = motif->get_residue_count();
	const int scaf_rescount = scaffold->get_residue_count();
	PRODART::POSE::int_vector res_mapping(mot_rescount, std::numeric_limits<int>::max());
	for (int m = 0; m < mot_rescount; m++){
		atom_shared_ptr m_ca_at = motif->get_bb_atom(POSE::CA, m);
		const vector3d m_ca_vec = m_ca_at->get_coords();
		double best_dist = std::numeric_limits<double>::max();
		for (int s = 0; s < scaf_rescount; s++){
			atom_shared_ptr s_ca_at = scaffold->get_bb_atom(POSE::CA, s);
			const vector3d s_ca_vec = s_ca_at->get_coords();
			const double dist = (m_ca_vec - s_ca_vec).mod();
			if (dist < best_dist){
				best_dist = dist;
				res_mapping[m] = s;
			}
		}
	}
	for (int m = 0; m < mot_rescount; m++){
		std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;
		if (res_mapping[m] < std::numeric_limits<int>::max()){
			residue_shared_ptr m_res = motif->get_residue(m);
			residue_shared_ptr s_res = scaffold->get_residue(res_mapping[m]);

			cout << "mapping motif res:\t" << m_res->get_trimmed_pdb_residue_index() << "\t"
					<< m_res->get_chain()->getChainID() << "\t"
					<< "scaffold res:\t" << s_res->get_trimmed_pdb_residue_index() << "\t"
					<< s_res->get_chain()->getChainID() << "\n";

			if (m_res->get_bb_atom(POSE::N)->isActiveAndSet() && s_res->get_bb_atom(POSE::N)->isActiveAndSet())
				atom_mapping[s_res->get_bb_atom(POSE::N)] = m_res->get_bb_atom(POSE::N);
			if (m_res->get_bb_atom(POSE::CA)->isActiveAndSet() && s_res->get_bb_atom(POSE::CA)->isActiveAndSet())
				atom_mapping[s_res->get_bb_atom(POSE::CA)] = m_res->get_bb_atom(POSE::CA);
			if (m_res->get_bb_atom(POSE::C)->isActiveAndSet() && s_res->get_bb_atom(POSE::C)->isActiveAndSet())
				atom_mapping[s_res->get_bb_atom(POSE::C)] = m_res->get_bb_atom(POSE::C);
			if (m_res->get_bb_atom(POSE::CB)->isActiveAndSet() && s_res->get_bb_atom(POSE::CB)->isActiveAndSet())
				atom_mapping[s_res->get_bb_atom(POSE::CB)] = m_res->get_bb_atom(POSE::CB);

			atom_shared_ptr_vector atoms_to_move = m_res->get_all_atoms();
			get_rmsd_superpose(atoms_to_move, atom_mapping);

			if (m_res->get_chain()->isPeptide()){
				s_res->set_type(m_res->get_type());
				s_res->get_sidechain()->copy_sidechain(m_res->get_sidechain());
				if (m_res->get_bb_atom(POSE::CB)->isActiveAndSet()
						&& !s_res->get_bb_atom(POSE::CB)->isActiveAndSet()){
					s_res->get_bb_atom(POSE::CB)->set_coords(m_res->get_bb_atom(POSE::CB)->get_coords());
					s_res->get_bb_atom(POSE::CB)->setActive(true);
					s_res->get_bb_atom(POSE::CB)->setSet(true);
				}
				scaffold->index();
			}

		}
		else {
			residue_shared_ptr m_res = motif->get_residue(m);
			cout << "no mapping found for motif res:\t" << m_res->get_trimmed_pdb_residue_index() << "\t"
					<< m_res->get_chain()->getChainID() << "\n";
		}
	}

	if (in_motif->get_all_atom_count() == motif->get_all_atom_count()){
		const int atom_count = in_motif->get_all_atom_count() ;
		std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;
		for ( int i = 0; i < atom_count; i++){
			atom_mapping[in_motif->get_atom(i)] =  motif->get_atom(i);
		}
		const double rmsd = get_rmsd(atom_mapping);
		cout << "final RMSD: " << rmsd << endl;
	}
	else {
		cerr << "lost atoms somewhere\n";
	}

}

PRODART::POSE::pose_shared_ptr make_centred_ALA_residue(){

	const sidechain_builder* sc_build = sidechain_builder::Instance();

	const_residue_reconstructor_shared_ptr res_constr =  sc_build->get_reconstructor(residue_type("ALA"));//new_residue_reconstructor();

	pose_shared_ptr blank_pose = new_pose();

	const double c_ca_cb_angle = UTILS::degrees_to_radians(111.0900);
	const double c_ca_bond = 1.5390;
	const double c_x = -c_ca_bond * sin(c_ca_cb_angle - (PI/2));
	const double c_y = c_ca_bond * cos(c_ca_cb_angle - (PI/2));

	blank_pose->add_new_chain('A', true);
	blank_pose->add_cterm_residue(residue_type("ALA"), 0);
	blank_pose->add_new_atom(vector3d(0,0,0), atom_type("CA"),0);
	blank_pose->add_new_atom(vector3d(1.5461,0,0), atom_type("CB"),0);
	blank_pose->add_new_atom(vector3d(c_x,c_y,0), atom_type("C"),0);
	blank_pose->index();

	res_constr->reconstruct_missing_atoms(blank_pose, 0);

	return blank_pose;
}

PRODART::POSE::pose_shared_ptr make_centred_trans_pept_bond(){
	pose_shared_ptr blank_pose = make_centred_ALA_residue();
	//const sidechain_builder* sc_build = sidechain_builder::Instance();
	//const_residue_reconstructor_shared_ptr res_constr =  sc_build->get_reconstructor(residue_type("ALA"));
	blank_pose->add_cterm_residue(residue_type("GLY"), 0);
	const vector3d n_pos = UTILS::dihedralEnd(blank_pose->get_bb_atom(POSE::N, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::CA, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			1.3558,
			degrees_to_radians(116.8400),
			degrees_to_radians(180.0000));
	const vector3d ca_pos = UTILS::dihedralEnd(blank_pose->get_bb_atom(POSE::CA, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			n_pos,
			1.4613,
			degrees_to_radians(126.7700),
			degrees_to_radians(180.0000));

	const vector3d o_pos = UTILS::dihedralEnd( n_pos,
			blank_pose->get_bb_atom(POSE::CA, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			1.2297,
			degrees_to_radians(122.5200),
			degrees_to_radians(180.0000));

	const vector3d h_pos = UTILS::dihedralEnd( blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			ca_pos,
			n_pos,
			1.2297,
			degrees_to_radians(122.5200),
			degrees_to_radians(180.0000));

	blank_pose->add_new_atom(o_pos, atom_type("O"),0);
	blank_pose->add_new_atom(n_pos, atom_type("N"),1);
	blank_pose->add_new_atom(ca_pos, atom_type("CA"),1);
	blank_pose->add_new_atom(h_pos, atom_type("H"),1);

	blank_pose->index();
	//res_constr->reconstruct_missing_atoms(blank_pose, 0);
	//res_constr->reconstruct_missing_atoms(blank_pose, 1);
	return blank_pose;
}
PRODART::POSE::pose_shared_ptr make_centred_cis_pept_bond(){
	pose_shared_ptr blank_pose = make_centred_ALA_residue();
	blank_pose->add_cterm_residue(residue_type("PRO"), 0);
	const vector3d n_pos = UTILS::dihedralEnd(blank_pose->get_bb_atom(POSE::N, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::CA, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			1.3569,
			degrees_to_radians(114.7500),
			degrees_to_radians(180.0000));
	const vector3d ca_pos = UTILS::dihedralEnd(blank_pose->get_bb_atom(POSE::CA, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			n_pos,
			1.4517,
			degrees_to_radians(124.8900),
			degrees_to_radians(0));

	const vector3d o_pos = UTILS::dihedralEnd( n_pos,
			blank_pose->get_bb_atom(POSE::CA, 0)->get_coords(),
			blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			1.2316,
			degrees_to_radians(120.4600),
			degrees_to_radians(177.1500));

	const vector3d h_pos = UTILS::dihedralEnd( blank_pose->get_bb_atom(POSE::C, 0)->get_coords(),
			ca_pos,
			n_pos,
			1.2297,
			degrees_to_radians(122.5200),
			degrees_to_radians(180.0000));

	blank_pose->add_new_atom(o_pos, atom_type("O"),0);
	blank_pose->add_new_atom(n_pos, atom_type("N"),1);
	blank_pose->add_new_atom(ca_pos, atom_type("CA"),1);
	blank_pose->add_new_atom(h_pos, atom_type("H"),1);

	blank_pose->index();
	return blank_pose;
}

PRODART::POSE::pose_shared_ptr copy_backbone_segment(PRODART::POSE::const_pose_shared_ptr protein, const int start, const int end,
		const bool renumber,
		const bool minimal_bb_set,
		const bool include_sc){
	pose_shared_ptr blank_pose = new_pose();
	if (start < 0 || start >= protein->get_residue_count() || end >= protein->get_residue_count() || end < start){
		cerr << "copy_backbone_segment: ERROR: start and end residues out of range" << endl;
		return blank_pose;
	}
	else if (protein->get_residue(start)->get_chain() != protein->get_residue(end)->get_chain()){
		cerr << "copy_backbone_segment: ERROR: start and end residues not on same chain" << endl;
		return blank_pose;
	}
	const char chainid = renumber ? 'A' : protein->get_residue(start)->get_chain()->getChainID();
	chain_shared_ptr nchain = blank_pose->add_new_chain(chainid, true);

	std::vector<BBAtomType> bb_types = minimal_bb_set ? prot_backbone_map::Instance()->get_smaller_bb_atom_list() : prot_backbone_map::Instance()->get_full_bb_atom_list();

	for (int i = start; i <= end; i++ ){
		const_residue_shared_ptr res = protein->get_residue(i);
		residue_shared_ptr nres = blank_pose->add_cterm_residue(res->get_type(), 0);
		if (!renumber) nres->set_pdb_residue_index(res->get_trimmed_pdb_residue_index());
		for (unsigned int at = 0; at < bb_types.size(); at++){
			const_atom_shared_ptr old_at = res->get_bb_atom(bb_types[at]);
			if (old_at->isActiveAndSet()){
				blank_pose->add_new_atom(old_at->get_coords(), old_at->get_type(),nres->get_internal_residue_index());
			}
		}
		if (include_sc){
			nres->set_sidechain(res->get_sidechain()->clone());
		}
	}



	blank_pose->index();


	return blank_pose;
}

bool output_centred_backbone_segments(PRODART::POSE::pose_shared_ptr protein, const int fraglen, std::ostream& output,
		const bool not_gly,
		const bool not_pro){

	if (!POSE_UTILS::validate_bb_pose(protein)){
		cerr << "output_centred_backbone_segments: WARNING: could not validate input pose:\t" << protein->get_label() << endl;
		return false;
	}
	else if (protein->get_residue_count() < fraglen){

	}

	const bool is_even = fraglen % 2 == 0 ? true : false;

	const int centre_res = fraglen / 2;

	pose_shared_ptr ideal_ala = PROTOCOLS::MISC::make_centred_ALA_residue();
	pose_shared_ptr ideal_trans = PROTOCOLS::MISC::make_centred_trans_pept_bond();
	pose_shared_ptr ideal_cis = PROTOCOLS::MISC::make_centred_cis_pept_bond();

	for (int i = 0; i<protein->get_residue_count()-(fraglen-1); i++){
		if (protein->get_residue(i)->get_chain() == protein->get_residue(i+(fraglen-1))->get_chain()){
			pose_shared_ptr bb_seg = PROTOCOLS::MISC::copy_backbone_segment(protein, i, i+(fraglen-1));
			if (POSE_UTILS::validate_bb_pose(bb_seg) && bb_seg->get_residue_count() == fraglen){
				const bool gly_skip = not_gly && bb_seg->get_residue(centre_res)->get_type() == residue_type("GLY");
				const bool pro_skip = not_pro && bb_seg->get_residue(centre_res)->get_type() == residue_type("PRO");
				if (!gly_skip && !pro_skip){
					if (!is_even){
						std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;
						atom_mapping[ideal_ala->get_bb_atom(POSE::N,0)] = bb_seg->get_bb_atom(POSE::N, centre_res);
						atom_mapping[ideal_ala->get_bb_atom(POSE::CA,0)] = bb_seg->get_bb_atom(POSE::CA, centre_res);
						atom_mapping[ideal_ala->get_bb_atom(POSE::C,0)] = bb_seg->get_bb_atom(POSE::C, centre_res);
						if (bb_seg->get_bb_atom(POSE::CB, centre_res)->isActiveAndSet()) atom_mapping[ideal_ala->get_bb_atom(POSE::CB,0)] = bb_seg->get_bb_atom(POSE::CB, centre_res);
						get_rmsd_superpose(bb_seg, atom_mapping);
						for (int j = 0; j <fraglen;j++){
							if (j != centre_res) bb_seg->get_bb_atom(POSE::CB, j)->setSet(false);
						}
					}
					else {
						// align peptide bond instead
						const double omega = bb_seg->get_omega_to_prev(centre_res);
						const bool is_cis = ( omega < (UTILS::PI / 2.0) && omega > -(UTILS::PI / 2.0)) ? true : false;
						pose_shared_ptr to_align = is_cis ? ideal_cis : ideal_trans ;
						std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;
						atom_mapping[to_align->get_bb_atom(POSE::CA,0)] = bb_seg->get_bb_atom(POSE::CA, centre_res-1);
						atom_mapping[to_align->get_bb_atom(POSE::C,0)] = bb_seg->get_bb_atom(POSE::C, centre_res-1);
						atom_mapping[to_align->get_bb_atom(POSE::O,0)] = bb_seg->get_bb_atom(POSE::O, centre_res-1);
						atom_mapping[to_align->get_bb_atom(POSE::N,1)] = bb_seg->get_bb_atom(POSE::N, centre_res);
						atom_mapping[to_align->get_bb_atom(POSE::CA,1)] = bb_seg->get_bb_atom(POSE::CA, centre_res);
						get_rmsd_superpose(bb_seg, atom_mapping);
						for (int j = 0; j <fraglen;j++){
							bb_seg->get_bb_atom(POSE::CB, j)->setSet(false);
						}
					}
					bb_seg->suppress_remark_output();
					bb_seg->outputPdb(output);
				}
			}
			else {
				cerr << "output_centred_backbone_segments: ERROR: could not validate backbone segment - probably a bug:\t" << protein->get_label() << endl;
				return false;
			}
		}
	}

	return true;
}
bool output_centred_backbone_segments(std::istream& filelist,
		const int fraglen,
		std::ostream& output,
		const bool not_gly,
		const bool not_pro){

	if (fraglen < 3){
		cerr << "output_centred_backbone_segments: ERROR: fraglen is less than 3: " << fraglen << endl;
		return false;
	}

	sidechain_builder::Instance();


    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;
    int pdbs_loaded = 0;

	cout << "reading file..." << endl;
    while ( !filelist.eof() ) {
            getline(filelist, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

                split( SplitVec, lineStr, is_any_of("\t ") );
                string this_pdb = SplitVec[0];

            	PRODART::POSE::pose_shared_ptr in_pose;

            	bool loadOK = true;
            	try {
            		cout << "output_centred_backbone_segments: parsing: " << this_pdb << " ... ";
            		in_pose = verified_load_bb_pose(this_pdb);//PRODART::POSE::new_pose();
            	}
            	catch (std::exception &e) {
        			cout << "failed" << endl;
            		std::cerr << "output_centred_backbone_segments: WARNING: skipping file: " << lineStr << endl;
            		loadOK = false;
            	}

            	if (loadOK){
            		if (output_centred_backbone_segments(in_pose, fraglen, output, not_gly, not_pro)){
            			cout << "done" << endl;
            			pdbs_loaded++;
            		}
            		else {
            			cout << "failed" << endl;
            		}
            	}



            }

    }

    cout << "output_centred_backbone_segments: INFO: " << pdbs_loaded << " PDB file has been successfully parsed" << endl;
    return true;
}


bool output_quality_filtered_list(std::istream& filelist,
		std::ostream& output){




    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;
    //int pdbs_loaded = 0;

	//cout << "reading file..." << endl;
    while ( !filelist.eof() ) {
            getline(filelist, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

                split( SplitVec, lineStr, is_any_of("\t ") );
                string this_pdb = SplitVec[0];

            	PRODART::POSE::pose_shared_ptr in_pose;

            	bool loadOK = true;
            	try {
            		//cout << "output_centred_backbone_segments: parsing: " << this_pdb << " ... ";
            		in_pose = verified_load_bb_pose(this_pdb);//PRODART::POSE::new_pose();
            	}
            	catch (std::exception &e) {
        			//cout << "failed" << endl;
            		std::cerr << "output_quality_filtered_list: INFO: PDB file failed: " << lineStr << endl;
            		loadOK = false;
            	}

            	if (loadOK){
            		output << lineStr << endl;
            	}



            }

    }

    return true;
}

PRODART::POSE::pose_shared_ptr get_single_chain_pose_byChainID(PRODART::POSE::const_pose_shared_ptr pose_, char chainID){
	const_chain_shared_ptr chain = pose_->get_chain(chainID);
	pose_shared_ptr sgl_chain_pose = new_pose();
	if (!chain){
		return sgl_chain_pose;
	}
	else {
		sgl_chain_pose->add_duplicated_chain(chain);
	}
	return sgl_chain_pose;

}

//! prefered function
PRODART::POSE::pose_shared_ptr get_single_chain_pose_byChainNum(PRODART::POSE::const_pose_shared_ptr pose_, const int chainnum){
	const_chain_shared_ptr chain = pose_->get_chain(chainnum);
	pose_shared_ptr sgl_chain_pose = new_pose();
	if (!chain){
		return sgl_chain_pose;
	}
	else {
		sgl_chain_pose->add_duplicated_chain(chain);
	}
	return sgl_chain_pose;
}

bool process_list_for_rosetta_match(std::istream& filelist,
		const std::string& output_dir,
		const int min_len,
		const int max_len){


	if (min_len > max_len){
		cerr << "process_pdb_file_for_rosetta_match: ERROR: min_len is more than max_len" << endl;
		return false;
	}

	path out_path(output_dir);

	if (exists(out_path)){
		if (is_directory(out_path)){

			string lineStr;
			unsigned long length, lineNum = 0 ;
			string_vector SplitVec;
			//int pdbs_loaded = 0;

			//cout << "reading file..." << endl;
			while ( !filelist.eof() ) {
				getline(filelist, lineStr);

				trim(lineStr);
				lineNum++;

				length = lineStr.length();

				//cout << endl << lineNum << " " << length << " ";

				if (length > 0) {

					split( SplitVec, lineStr, is_any_of("\t ") );
					string this_pdb = SplitVec[0];

					process_pdb_file_for_rosetta_match(this_pdb, output_dir, min_len, max_len);



				}

			}

			return true;
		}
		else {
			cerr << "process_pdb_file_for_rosetta_match: ERROR: output directory is not a directory" << endl;
			return false;
		}
	}
	else {
		cerr << "process_pdb_file_for_rosetta_match: ERROR: output directory does not exist" << endl;
		return false;
	}

	return false;
}

bool process_pdb_file_for_rosetta_match(std::string& pose_file,
		const std::string& output_dir,
		const int min_len,
		const int max_len){

	if (min_len > max_len){
		cerr << "process_pdb_file_for_rosetta_match: ERROR: min_len is more than max_len" << endl;
		return false;
	}

	PRODART::POSE::pose_shared_ptr in_pose;

	bool loadOK = true;
	try {
		//cout << "output_centred_backbone_segments: parsing: " << this_pdb << " ... ";
		in_pose = verified_load_bb_pose(pose_file);//PRODART::POSE::new_pose();
	}
	catch (std::exception &e) {
		//cout << "failed" << endl;
		std::cerr << "process_pdb_file_for_rosetta_match: INFO: PDB file failed: " << pose_file << endl;
		loadOK = false;
	}

	if (loadOK){
		// do stuff here
		cout << pose_file << endl;


		path out_path(output_dir);
		path input_flpth(pose_file);
		path out_stem = out_path / input_flpth.stem();
		if (exists(out_path)){
			if (is_directory(out_path)){
				const int chain_count = in_pose->get_chain_count();
				for (int i = 0; i< chain_count; i++){
					chain_shared_ptr this_chain = in_pose->get_chain(i);
					if (this_chain->isPeptide()){
						const char chainid = this_chain->getChainID();
						pose_shared_ptr sgl_chain_pose = new_pose();
						sgl_chain_pose->add_duplicated_chain(this_chain);
						if (validate_bb_pose(sgl_chain_pose)){
							if (sgl_chain_pose->get_residue_count() <= max_len
									&& sgl_chain_pose->get_residue_count() >= min_len ){
								string out_tail(out_stem.string());
								out_tail.append("_");
								out_tail.append(string(1,chainid));
								out_tail.append(".pdb");
								path out_file(out_tail);
								cout << "outputting file : " << out_file.string() << endl;
								boost::filesystem::ofstream output(out_file);
								if (output.is_open()){
									sgl_chain_pose->outputPdb(output);
								}
							}
							else {
								// too short
								std::cerr << "process_pdb_file_for_rosetta_match: INFO: PDB chain: " << chainid
										<< " is too short or too long: " << pose_file
										<< " length: " <<  sgl_chain_pose->get_residue_count()
										<< endl;
							}
						}
						else {
							std::cerr << "process_pdb_file_for_rosetta_match: ERROR: Probably bug: PDB chain: " << chainid << " could not be validated: " << pose_file << endl;

						}

					}

				}
			}
			else {
				cerr << "process_pdb_file_for_rosetta_match: ERROR: output directory is not a directory" << endl;
				return false;
			}
		}
		else {
			cerr << "process_pdb_file_for_rosetta_match: ERROR: output directory does not exist" << endl;
			return false;
		}
	}
	else {
		return false;
	}

	return true;
}


bool output_design_quality(std::istream& filelist,
		 std::ostream& output){

	string lineStr;
	unsigned long length, lineNum = 0 ;
	string_vector SplitVec;

	while ( !filelist.eof() ) {
		getline(filelist, lineStr);

		trim(lineStr);
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {

			split( SplitVec, lineStr, is_any_of("\t ") );
			if (SplitVec.size() >= 2){
				string des_pdb = SplitVec[0];
				string orig_scaff = SplitVec[1];

				cout << "des: " << des_pdb << " scaff: " << orig_scaff << endl;


				pose_shared_ptr des_pose, scaff_pose;

				bool load_ok = true;

				try{
					des_pose = unverified_load_pose(des_pdb);

				}
				catch (std::exception &e) {
					//cout << "failed" << endl;
					std::cerr << "output_design_quality: INFO: PDB file failed: " << des_pdb << endl;
					load_ok = false;
				}

				try {
					scaff_pose = unverified_load_pose(orig_scaff);
				}
				catch (std::exception &e) {
					//cout << "failed" << endl;
					std::cerr << "output_design_quality: INFO: PDB file failed: " << orig_scaff << endl;
					load_ok = false;
				}

				if (load_ok
						&& des_pose
						&& scaff_pose){

					const bool des_ok = validate_bb_pose(des_pose);
					const bool scaff_ok = validate_bb_pose(scaff_pose);

					const bool des_high_tol_ok = validate_bb_pose(des_pose, 0.15);

					// do these pdbs match
					bool pdbs_match = true;

					if (des_pose->get_residue_count() != scaff_pose->get_residue_count()){
						pdbs_match = false;
					}
					if (des_pose->get_chain_count() != scaff_pose->get_chain_count()){
						pdbs_match = false;
					}
					else {
						for (int i = 0; i < des_pose->get_chain_count(); i++){
							if (des_pose->get_chain(i)->length() != scaff_pose->get_chain(i)->length()){
								pdbs_match = false;
							}
						}
					}

					if (pdbs_match){
						// carry on with scoring


						boost::tuple<double, double, double> max_moves = get_max_alpha_beta_loop_N_CA_C_movement(scaff_pose, des_pose );
						const int non_cisPro_conversions = get_cisPro_to_nonPro_count(scaff_pose, des_pose );
						const bool des_worse_than_scaff = (scaff_ok == true) && (des_ok == false);

						output << des_pdb << "\t" << orig_scaff << "\t" << max_moves.get<0>() << "\t"
								<< max_moves.get<1>() << "\t"
								<< max_moves.get<2>() << "\t"
								<< des_worse_than_scaff << "\t"
								<< !des_high_tol_ok << "\t"
								<< non_cisPro_conversions << "\t"
								<< endl;

					}
					else{
						cout << des_pdb << "\t" << orig_scaff << "\t" << "match_error:skipping" << endl;
					}

				}
				else {
					cout << des_pdb << "\t" << orig_scaff << "\t" << "load_error:skipping" << endl;
				}
			}

		}
	}
	return true;
}


boost::tuple<double, double, double> get_max_alpha_beta_loop_N_CA_C_movement(const PRODART::POSE::pose_shared_ptr ref,
		const PRODART::POSE::pose_shared_ptr protein){


	if (ref->get_residue_count() != protein->get_residue_count()){
		return boost::make_tuple(double(0.0),double(0.0),double(0.0));//tuple<double, double, double>(0.0,0.0,0.0) ;
	}

	POSE_UTILS::quick_add_HN(ref,false);
	for (int i = 0; i < ref->get_residue_count(); i++){
		atom_shared_ptr at = ref->get_bb_atom(POSE::H, i);
		at->set_type(atom_type("H"));
	}
	POSE::META::bb_pose_meta_shared_ptr meta = new_bb_pose_meta(ref);
	const three_state_sec_struct_vector secs = meta->get_sec_struct();
	cout << secs << endl;

	double max_alpha = 0, max_beta = 0, max_loop = 0;


	for (int i = 0; i < ref->get_residue_count(); i++){
		const double n_dist = (ref->get_bb_atom(POSE::N, i)->get_coords() - protein->get_bb_atom(POSE::N, i)->get_coords()).mod();
		const double ca_dist = (ref->get_bb_atom(POSE::CA, i)->get_coords() - protein->get_bb_atom(POSE::CA, i)->get_coords()).mod();
		const double c_dist = (ref->get_bb_atom(POSE::CA, i)->get_coords() - protein->get_bb_atom(POSE::CA, i)->get_coords()).mod();

		const double max_val = max(n_dist, max(ca_dist, c_dist));

		if (secs[i] == ss3_HELIX){
			max_alpha = max(max_alpha,max_val );
		}
		else if (secs[i] == ss3_STRAND){
			max_beta = max(max_beta, max_val);
		}
		else {
			max_loop = max(max_loop, max_val);
		}


	}


	return boost::make_tuple(max_alpha, max_beta, max_loop);
}


int get_cisPro_to_nonPro_count(const PRODART::POSE::pose_shared_ptr ref,
		const PRODART::POSE::pose_shared_ptr protein){
	if (ref->get_residue_count() != protein->get_residue_count()){
		return -1;
	}

	const four_state_sec_struct_vector ref_secs = POSE_UTILS::get_phi_psi_omega_sector(ref);
	const four_state_sec_struct_vector prot_secs = POSE_UTILS::get_phi_psi_omega_sector(protein);

	int count = 0;

	for (int i = 0; i < ref->get_residue_count(); i++){
		if (ref_secs[i] == ss4_CIS && ref->get_residue(i)->get_type().is_equal3(residue_type("PRO"))){
			if (!protein->get_residue(i)->get_type().is_equal3(residue_type("PRO"))){
				count++;
			}
		}
	}


	return count;
}


bool output_seq_diff( std::ostream& output,
		PRODART::POSE::pose_shared_ptr ref,
		PRODART::POSE::pose_shared_ptr protein){
	if (ref->get_residue_count() != protein->get_residue_count()){
		return false;
	}
	PRODART::POSE::int_vector diff_res;
	for (int i = 0; i < ref->get_residue_count(); i++){
		if (!ref->get_residue(i)->get_type().is_equal3(
				protein->get_residue(i)->get_type())){
			diff_res.push_back(i);
		}
	}

	output << "resid ";
	for (unsigned int i = 0; i < diff_res.size(); i++){
		output << protein->get_residue(diff_res[i])->get_trimmed_pdb_residue_index() << " ";
	}
	output << endl << endl;

	output << "resid ";
	for (unsigned int i = 0; i < diff_res.size(); i++){
		output << protein->get_residue(diff_res[i])->get_trimmed_pdb_residue_index() << "+";
	}
	output << endl << endl;

	for (unsigned int i = 0; i < diff_res.size(); i++){
		output << ref->get_residue(diff_res[i])->get_type().get_label3() << "\t"
				<< protein->get_residue(diff_res[i])->get_trimmed_pdb_residue_index() << "\t"
				<< protein->get_residue(diff_res[i])->get_type().get_label3() << "\t"
				<< endl;
	}


	return true;
}


bool output_seq3( std::ostream& output,
		PRODART::POSE::pose_shared_ptr protein){
	for (int i = 0; i < protein->get_residue_count(); i++){
		output	<< protein->get_residue(i)->get_trimmed_pdb_residue_index() << "\t"
				<< protein->get_residue(i)->get_chain()->getChainID() << "\t"
				<< protein->get_residue(i)->get_type().get_label3() << "\t"
				<< endl;
	}

	return true;

}

bool output_seq1( std::ostream& output,
		PRODART::POSE::pose_shared_ptr protein){
	const int chain_count = protein->get_chain_count();
	for (int ch = 0 ; ch < chain_count; ch++){
		chain_shared_ptr chain_ = protein->get_chain(ch);
		residue_type_vector ch_seq = chain_->get_sequence();
		string seqStr =  POSE::residue_type_map::Instance()->to_1_letter_code(ch_seq);
		output  << "> " << protein->get_label() << "|" << chain_->getChainID() << "\n"
				<< seqStr
				<< endl;
	}
	return true;
}


PRODART::POSE::atom_shared_ptr_vector get_atom_subset(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask, PRODART::POSE::atom_type_vector &vec){
	atom_shared_ptr_vector moved_atoms;


	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){


			for (unsigned int j = 0; j < vec.size(); j++){
				if (protein->get_atom(vec[j], i)->isActive()) moved_atoms.push_back(protein->get_atom(vec[j], i));
			}

		}
	}

	return moved_atoms;
}


PRODART::POSE::atom_shared_ptr_vector get_atom_subset(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask, std::vector<PRODART::POSE::BBAtomType> &vec){
	atom_shared_ptr_vector moved_atoms;


	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){


			for (unsigned int j = 0; j < vec.size(); j++){
				if (protein->get_bb_atom(vec[j], i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(vec[j], i));
			}

		}
	}

	return moved_atoms;
}


void output_raw_bond_lengths_angles(PRODART::POSE::const_pose_shared_ptr protein, std::ostream& output){


	const int rescount = protein->get_residue_count();

	for (int i = 0 ; i < rescount; i++){

		const_residue_shared_ptr res  = protein->get_residue(i);
		const_chain_shared_ptr ch = res->get_chain();

		if (ch->isPeptide()){

			const_atom_shared_ptr c_m1 = const_atom_shared_ptr();
			if (res->get_prev_residue()) c_m1 = res->get_prev_residue()->get_bb_atom(C);
			const_atom_shared_ptr n = protein->get_bb_atom(N, i);
			const_atom_shared_ptr ca = protein->get_bb_atom(POSE::CA, i);
			const_atom_shared_ptr c = protein->get_bb_atom(POSE::C, i);
			const_atom_shared_ptr o = protein->get_bb_atom(POSE::O, i);
			const_atom_shared_ptr n_p1 = const_atom_shared_ptr();
			if (res->get_next_residue()) n_p1 = res->get_next_residue()->get_bb_atom(N);

			const double n_ca_c = angle(n->get_coords(), ca->get_coords(), c->get_coords());
			const double c_n_ca = c_m1 ? angle(c_m1->get_coords(), n->get_coords(), ca->get_coords()) : degrees_to_radians(999.999);
			const double ca_c_n = n_p1 ? angle(ca->get_coords(), c->get_coords(), n_p1->get_coords()) : degrees_to_radians(999.999);
			const double n_c_o = n_p1 ? angle(n_p1->get_coords(), c->get_coords(), o->get_coords()) : degrees_to_radians(999.999);

			const double cm1_n = c_m1 ? c_m1->get_coords().dist(n->get_coords()) : 0.0;
			const double n_ca = n->get_coords().dist(ca->get_coords()) ;
			const double ca_c = ca->get_coords().dist(c->get_coords()) ;
			const double c_o = o->get_coords().dist(c->get_coords()) ;
			//const double c_np1 =  n_p1 ? c->get_coords().dist(n_p1->get_coords()) : 0.0;

			const char chid = (res->get_chain()->getChainID() != ' ') ? res->get_chain()->getChainID() : '-';

			output  << fixed << setprecision(3)
					<< protein->get_label() << "\t"
					<< res->get_trimmed_pdb_residue_index() << "\t"
					<< chid << "\t"
					<< res->get_type().get_label3() << "\t"
					<< radians_to_degrees(c_n_ca) << "\t"
					<< radians_to_degrees(n_ca_c) << "\t"
					<< radians_to_degrees(ca_c_n) << "\t"
					<< radians_to_degrees(n_c_o) << "\t"
					<< setprecision(4)
					<< cm1_n << "\t"
					<< n_ca << "\t"
					<< ca_c << "\t"
					//<< c_np1 << "\t"
					<< c_o << "\t"
					<< (!res->is_terminal() ? "middle" : "terminal")
					<< endl;

		}

	}

}


void output_raw_bond_lengths_angles( std::istream& filelist,
		std::ostream& output){


    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;
    //int pdbs_loaded = 0;

	//cout << "reading file..." << endl;
    while ( !filelist.eof() ) {
            getline(filelist, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

                split( SplitVec, lineStr, is_any_of("\t ") );
                string this_pdb = SplitVec[0];

            	PRODART::POSE::pose_shared_ptr in_pose = new_pose();

            	std::ifstream protein_file(this_pdb.c_str(), ios::in);

            	bool loadOK = true;

            	if (protein_file.is_open()){
            		in_pose->loadPdb(protein_file);
            		in_pose->set_label(this_pdb);
            	}
            	else {
            		loadOK = false;
            		std::cerr << "output_raw_bond_lengths_angles: INFO: PDB file failed to load: " << lineStr << endl;

            	}


            	/*
            	try {
            		in_pose = verified_load_bb_pose(this_pdb);//PRODART::POSE::new_pose();
            	}
            	catch (std::exception &e) {
            		std::cerr << "output_quality_filtered_list: INFO: PDB file failed: " << lineStr << endl;
            		loadOK = false;
            	}
            	*/

            	if (loadOK){
            		output_raw_bond_lengths_angles(in_pose, output);
            	}



            }

    }
}


double get_max_motif_atom_deviation(const int start_res, PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		double& rtn_rmsd){
	pose_shared_ptr motif_cp = motif->clone();

	bool isError = false;
	std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;
	for (int i = 0; i < motif_cp->get_residue_count(); i++ ){
		if (i + start_res < protein->get_residue_count()){
			const_residue_shared_ptr  prot_res = protein->get_residue(i + start_res);
			const_residue_shared_ptr mot_res = motif_cp->get_residue(i);

			if(prot_res && mot_res){

				{
					const_atom_shared_ptr N_mot = mot_res->get_bb_atom(POSE::N);
					const_atom_shared_ptr N_prot = prot_res->get_bb_atom(POSE::N);
					if (N_mot->isActiveAndSet() && N_prot->isActiveAndSet()){
						atom_mapping[N_prot] = N_mot;
					}
					else {
						isError = true;
					}
				}

				{
					const_atom_shared_ptr CA_mot = mot_res->get_bb_atom(POSE::CA);
					const_atom_shared_ptr CA_prot = prot_res->get_bb_atom(POSE::CA);
					if (CA_mot->isActiveAndSet() && CA_prot->isActiveAndSet()){
						atom_mapping[CA_prot] = CA_mot;
					}
					else {
						isError = true;
					}
				}

				{
					const_atom_shared_ptr C_mot = mot_res->get_bb_atom(POSE::C);
					const_atom_shared_ptr C_prot = prot_res->get_bb_atom(POSE::C);
					if (C_mot->isActiveAndSet() && C_prot->isActiveAndSet()){
						atom_mapping[C_prot] = C_mot;
					}
					else {
						isError = true;
					}
				}

				{
					const_atom_shared_ptr O_mot = mot_res->get_bb_atom(POSE::O);
					const_atom_shared_ptr O_prot = prot_res->get_bb_atom(POSE::O);
					if (O_mot->isActiveAndSet() && O_prot->isActiveAndSet()){
						atom_mapping[O_prot] = O_mot;
					}
					else {
						isError = true;
					}
				}

			}
			else {
				isError = true;
				break;
			}
		}
		else {
			//run over
			isError = true;
			break;
		}
	}

	if (isError){
		return std::numeric_limits<double>::max();
	}


	const double rmsd = get_rmsd_superpose(motif_cp, atom_mapping);
	rtn_rmsd = rmsd;

	double max = 0;

	for (std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr>::iterator iter = atom_mapping.begin();
			iter != atom_mapping.end(); iter++){

		const double dev = (iter->first->get_coords() - iter->second->get_coords()).mod();
		if (dev > max){
			max = dev;
		}

	}
	return max;
}


void scan_motif_atom_deviation(POSE::double_vector& min_vals,
		std::vector<std::string>& mot_res,
		PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const bool use_rmsd ){
	//double_vector min_vals(protein->get_residue_count(), std::numeric_limits<double>::max());
	PRODART::POSE::int_vector min_res(protein->get_residue_count(), -1);
	double rmsd = std::numeric_limits<double>::max();

	for (int i = 0 ; i < protein->get_residue_count(); i++){
		double max_val = get_max_motif_atom_deviation(i, protein, motif, rmsd);
		if (use_rmsd){
			max_val = rmsd;
		}

		for (int j = i; j < i +motif->get_residue_count(); j++){
			if (j < protein->get_residue_count()){
				if (max_val < min_vals[j]){
					min_vals[j] = max_val;
					min_res[j] = j -i;

					const_residue_shared_ptr m_res;
					if (min_res[j] >= 0 && min_res[j] < motif->get_residue_count()){
						m_res = motif->get_residue(min_res[j]);
					}
					string resid = "UNK"; //

					if (m_res){
						resid = m_res->get_pdb_residue_index();
						trim(resid);
					}

					mot_res[j] = resid;
				}
			}
		}
	}

	/*
	for (int i = 0 ; i < protein->get_residue_count(); i++){

		const_residue_shared_ptr res = protein->get_residue(i);
		//const char chid = (res->get_chain()->getChainID() != ' ') ? res->get_chain()->getChainID() : '-';

	}
	*/
}

void scan_motif_atom_deviation(std::ostream& output,
		PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const bool use_rmsd ){

	PRODART::POSE::double_vector min_vals(protein->get_residue_count(), std::numeric_limits<double>::max());
	std::vector<std::string> mot_res(protein->get_residue_count(), "UNK");


	for (int i = 0; i < motif->get_chain_count(); i++){
		pose_shared_ptr ch_motif = get_single_chain_pose_byChainNum(motif, i);
		scan_motif_atom_deviation(min_vals, mot_res, protein, ch_motif, use_rmsd);
	}

	for (int i = 0 ; i < protein->get_residue_count(); i++){

		const_residue_shared_ptr res = protein->get_residue(i);
		const char chid = (res->get_chain()->getChainID() != ' ') ? res->get_chain()->getChainID() : '-';



		output  << fixed << setprecision(4)
				<< protein->get_label() << "\t"
				<< res->get_trimmed_pdb_residue_index() << "\t"
				<< chid << "\t"
				<< res->get_type().get_label3() << "\t"
				<< min_vals[i] << "\t"
				<< mot_res[i] << "\t"
				<< endl;
	}

}

void scan_motif_atom_deviation(std::ostream& output,
		 std::istream& filelist,
		PRODART::POSE::const_pose_shared_ptr motif,
		const bool use_rmsd ){

    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;
    //int pdbs_loaded = 0;

	//cout << "reading file..." << endl;
    while ( !filelist.eof() ) {
            getline(filelist, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

                split( SplitVec, lineStr, is_any_of("\t ") );
                string this_pdb = SplitVec[0];

            	PRODART::POSE::pose_shared_ptr in_pose = new_pose();

            	std::ifstream protein_file(this_pdb.c_str(), ios::in);

            	bool loadOK = true;

            	if (protein_file.is_open()){
            		in_pose->loadPdb(protein_file);
            		in_pose->set_label(this_pdb);
            	}
            	else {
            		loadOK = false;
            		std::cerr << "scan_motif_atom_deviation: INFO: PDB file failed to load: " << lineStr << endl;

            	}


            	/*
            	try {
            		in_pose = verified_load_bb_pose(this_pdb);//PRODART::POSE::new_pose();
            	}
            	catch (std::exception &e) {
            		std::cerr << "output_quality_filtered_list: INFO: PDB file failed: " << lineStr << endl;
            		loadOK = false;
            	}
            	*/

            	if (loadOK){
            		scan_motif_atom_deviation(output, in_pose, motif, use_rmsd);
            	}



            }

    }
}

bool is_gap_next(PRODART::POSE::const_pose_shared_ptr protein,
		const int resnum,
		const double cutoff){


	if (resnum < protein->get_residue_count()){

		const_residue_shared_ptr res = protein->get_residue(resnum);

		if (!res->get_next_residue()){
			return false;
		}

		const_atom_shared_ptr this_ca = protein->get_bb_atom(POSE::CA, resnum);
		const_atom_shared_ptr next_ca = protein->get_bb_atom(POSE::CA, resnum+1);
		if (this_ca && next_ca){
			if (this_ca->isSet() && next_ca->isSet()){
				const double dist = (this_ca->get_coords() - next_ca->get_coords()).mod();
				if (dist > cutoff){
					return true;
				}
			}
		}
	}
	else {
		return false;
	}

	return false;
}

//! outputs the lowest maximum backbone atom deviation in angstroms from the motif by sliding motif along the protein
void scan_motif_atom_deviation_annot_seq(std::ostream& output,
		PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd){

	PRODART::POSE::double_vector min_vals(protein->get_residue_count(), std::numeric_limits<double>::max());
	std::vector<std::string> mot_res(protein->get_residue_count(), "UNK");


	for (int i = 0; i < motif->get_chain_count(); i++){
		pose_shared_ptr ch_motif = get_single_chain_pose_byChainNum(motif, i);
		scan_motif_atom_deviation(min_vals, mot_res, protein, ch_motif, use_rmsd);
	}


	for (int ch_num = 0; ch_num < protein->get_chain_count(); ch_num++){

		const_chain_shared_ptr ch = protein->get_chain(ch_num);
		const char chid = (ch->getChainID() != ' ') ? ch->getChainID() : '-';

		output << ">" << protein->get_label() << " chain:" << chid
				<< endl;

		for (int i = ch->get_first_internal_residue_index() ; i <= ch->get_last_internal_residue_index(); i++){
			const_residue_shared_ptr res = protein->get_residue(i);
			output  << POSE::residue_type_map::Instance()->to_1_letter_code(res->get_type());

			if (is_gap_next(protein,i, 4.1)){
				output << "-";
			}
		}
		output << endl;

		for (int i = ch->get_first_internal_residue_index() ; i <= ch->get_last_internal_residue_index(); i++){
			const_residue_shared_ptr res = protein->get_residue(i);

			if (min_vals[i] < cutoff){
				output << mot_res[i];
			}
			else {
				output << " ";
			}

			if (is_gap_next(protein,i, 4.1)){
				output << " ";
			}

		}
		output << endl;


	}

}


void scan_motif_atom_deviation_annot_seq(std::ostream& output,
		 std::istream& filelist,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd){
    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;
    //int pdbs_loaded = 0;

	//cout << "reading file..." << endl;
    while ( !filelist.eof() ) {
            getline(filelist, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

                split( SplitVec, lineStr, is_any_of("\t ") );
                string this_pdb = SplitVec[0];

            	PRODART::POSE::pose_shared_ptr in_pose = new_pose();

            	std::ifstream protein_file(this_pdb.c_str(), ios::in);

            	bool loadOK = true;

            	if (protein_file.is_open()){
            		in_pose->loadPdb(protein_file);
            		in_pose->set_label(this_pdb);
            	}
            	else {
            		loadOK = false;
            		std::cerr << "scan_motif_atom_deviation: INFO: PDB file failed to load: " << lineStr << endl;

            	}


            	/*
            	try {
            		in_pose = verified_load_bb_pose(this_pdb);//PRODART::POSE::new_pose();
            	}
            	catch (std::exception &e) {
            		std::cerr << "output_quality_filtered_list: INFO: PDB file failed: " << lineStr << endl;
            		loadOK = false;
            	}
            	*/

            	if (loadOK){
            		scan_motif_atom_deviation_annot_seq(output, in_pose, motif, cutoff, use_rmsd);
            	}



            }

    }
}


//! outputs all fragments
void scan_motif_atom_deviation_output_frags(PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd ){

	PRODART::POSE::double_vector min_vals(protein->get_residue_count(), std::numeric_limits<double>::max());
	std::vector<std::string> mot_res(protein->get_residue_count(), "UNK");

	std::set<int> uniq_mot_nres;
	for (int i = 0 ; i < motif->get_residue_count(); i++){
		string rid = motif->get_residue(i)->get_pdb_residue_index();
		trim(rid);
		uniq_mot_nres.insert(lexical_cast<int>(rid));
	}
	const int mot_len = uniq_mot_nres.size();
	const int mot_multiples = 2;


	for (int i = 0; i < motif->get_chain_count(); i++){
		pose_shared_ptr ch_motif = get_single_chain_pose_byChainNum(motif, i);
		scan_motif_atom_deviation(min_vals, mot_res, protein, ch_motif, use_rmsd);
	}




	for (int ch_num = 0; ch_num < protein->get_chain_count(); ch_num++){

		const_chain_shared_ptr ch = protein->get_chain(ch_num);
		//const char chid = (ch->getChainID() != ' ') ? ch->getChainID() : '-';


		for (int i = ch->get_first_internal_residue_index() ; i <= ch->get_last_internal_residue_index(); i++){

			bool is_frag = true;
			for (int j = i; j < i+(mot_len*mot_multiples); j++){
				if (j >= ch->get_last_internal_residue_index()){
					is_frag = false;
					break;
				}

				//const_residue_shared_ptr res = protein->get_residue(j);
				if (min_vals[j] >= cutoff){
					is_frag = false;
					break;
				}


				if (is_gap_next(protein,j, 4.1)){
					is_frag = false;
					break;
				}
			}

			if (is_frag){
				pose_shared_ptr frag = copy_backbone_segment( protein,
						i,
						i+(mot_len*mot_multiples)-1,
						true,
						true,
						false);

				for (int rn = 0; rn < frag->get_residue_count(); rn++ ){
					frag->get_bb_atom(CB,rn)->setSet(false);
				}


				pose_shared_ptr ch_motif = new_pose();
				for (int chn = 0; chn < motif->get_chain_count(); chn++){
					if (mot_res[i].compare(motif->get_chain(chn)->get_residue(0)->get_trimmed_pdb_residue_index()) == 0){
						ch_motif = get_single_chain_pose_byChainNum(motif, chn);
						break;
					}
				}

				bool isError = false;
				std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;
				for (int rn = 0; rn < ch_motif->get_residue_count(); rn++ ){
					if (rn < frag->get_residue_count()){
						const_residue_shared_ptr  prot_res = frag->get_residue(rn);
						const_residue_shared_ptr mot_res = ch_motif->get_residue(rn);

						if(prot_res && mot_res){

							{
								const_atom_shared_ptr N_mot = mot_res->get_bb_atom(POSE::N);
								const_atom_shared_ptr N_prot = prot_res->get_bb_atom(POSE::N);
								if (N_mot->isActiveAndSet() && N_prot->isActiveAndSet()){
									atom_mapping[N_mot] = N_prot;
								}
								else {
									isError = true;
								}
							}

							{
								const_atom_shared_ptr CA_mot = mot_res->get_bb_atom(POSE::CA);
								const_atom_shared_ptr CA_prot = prot_res->get_bb_atom(POSE::CA);
								if (CA_mot->isActiveAndSet() && CA_prot->isActiveAndSet()){
									atom_mapping[CA_mot] = CA_prot;
								}
								else {
									isError = true;
								}
							}

							{
								const_atom_shared_ptr C_mot = mot_res->get_bb_atom(POSE::C);
								const_atom_shared_ptr C_prot = prot_res->get_bb_atom(POSE::C);
								if (C_mot->isActiveAndSet() && C_prot->isActiveAndSet()){
									atom_mapping[C_mot] = C_prot;
								}
								else {
									isError = true;
								}
							}

							{
								const_atom_shared_ptr O_mot = mot_res->get_bb_atom(POSE::O);
								const_atom_shared_ptr O_prot = prot_res->get_bb_atom(POSE::O);
								if (O_mot->isActiveAndSet() && O_prot->isActiveAndSet()){
									atom_mapping[O_mot] = O_prot;
								}
								else {
									isError = true;
								}
							}

						}
						else {
							isError = true;
							break;
						}
					}
					else {
						//run over
						isError = true;
						break;
					}
				}

				const double rmsd = isError==false ? get_rmsd_superpose(frag, atom_mapping) : 9999;

				cout << rmsd << endl;

				string outname("scan_motif_frag_");
				outname.append(mot_res[i]);
				outname.append(".pdb");


				POSE_UTILS::to_homopolymer(frag, residue_type("GLY"));
				ofstream output_pdb_file(outname.c_str(), ios::out | ios::app);
				frag->outputPdb(output_pdb_file);

			}

		}


	}

}


void scan_motif_atom_deviation_output_frags(std::istream& filelist,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd ){
    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;
    //int pdbs_loaded = 0;

	//cout << "reading file..." << endl;
    while ( !filelist.eof() ) {
            getline(filelist, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

                split( SplitVec, lineStr, is_any_of("\t ") );
                string this_pdb = SplitVec[0];

            	PRODART::POSE::pose_shared_ptr in_pose = new_pose();

            	std::ifstream protein_file(this_pdb.c_str(), ios::in);

            	bool loadOK = true;

            	if (protein_file.is_open()){
            		in_pose->loadPdb(protein_file);
            		in_pose->set_label(this_pdb);
            	}
            	else {
            		loadOK = false;
            		std::cerr << "scan_motif_atom_deviation: INFO: PDB file failed to load: " << lineStr << endl;

            	}


            	/*
            	try {
            		in_pose = verified_load_bb_pose(this_pdb);//PRODART::POSE::new_pose();
            	}
            	catch (std::exception &e) {
            		std::cerr << "output_quality_filtered_list: INFO: PDB file failed: " << lineStr << endl;
            		loadOK = false;
            	}
            	*/

            	if (loadOK){
            		scan_motif_atom_deviation_output_frags( in_pose, motif, cutoff, use_rmsd);
            	}



            }

    }


}



PRODART::POSE::pose_shared_ptr invert_chain_direction(PRODART::POSE::const_pose_shared_ptr protein){

	PRODART::POSE::pose_shared_ptr n_pose = POSE::new_pose();

	for (int ch = 0; ch < protein->get_chain_count(); ch++){
		const_chain_shared_ptr o_chain = protein->get_chain(ch);
		chain_shared_ptr n_chain = n_pose->add_new_chain(o_chain->getChainID(), true);
		for (int r = o_chain->length()-1; r >= 0; r--){
			residue_shared_ptr n_res = n_pose->add_cterm_residue(residue_type("ALA"), ch);
			n_pose->add_new_atom(o_chain->get_ca_pos(r),atom_type("CA"), n_res->get_internal_residue_index());
		}
		n_pose->index();
	}
	CA2MAIN::new_ca2main_minimise_fixed_ca(n_pose);
	return n_pose;
}

PRODART::POSE::pose_shared_ptr make_extended_ca_chain(const PRODART::POSE::residue_type_vector seq,
		const char chainID,
		const int start_res,
		const double random_component){

	MTRand::MTRand_shared_ptr randomNumGen = PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id());
	PRODART::POSE::pose_shared_ptr n_pose = POSE::new_pose();
	//const char ch = 'A';
	chain_shared_ptr n_chain = n_pose->add_new_chain(chainID, true);
	vector3d pos(0,0,0);
	vector3d incr(3.8,0,0);
	for (unsigned int r = 0; r < seq.size(); r++){
		residue_shared_ptr n_res = n_pose->add_cterm_residue(seq[r], 0);
		const vector3d randVec(randomNumGen->randNorm(0, random_component) ,
				randomNumGen->randNorm(0, random_component) ,
				randomNumGen->randNorm(0, random_component));
		n_pose->add_new_atom(pos+randVec,atom_type("CA"), r);
		pos = pos + incr;
	}
	n_pose->index();
    n_pose->renumber_residues(chainID, start_res);
	vector3d com = get_ca_CoM(n_pose);
	for (unsigned int r = 0; r < seq.size(); r++){
		n_pose->get_bb_atom(POSE::CA, r)->set_coords(n_pose->get_bb_atom(POSE::CA, r)->get_coords() - com );
	}

	bool_vector atom_selection(n_pose->get_all_atom_count(), true);
	LOOP::ca_minimise_bonds_bumps_angles(n_pose, atom_selection, false);


	return n_pose;
}

void extend_ca_chain(PRODART::POSE::pose_shared_ptr  protein){
	MTRand::MTRand_shared_ptr randomNumGen = PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id());
	const double random_component = 1.0;
	vector3d pos(0,0,0);
	vector3d incr(3.8,0,0);
	for ( int r = 0; r < protein->get_residue_count(); r++){
		const vector3d randVec(randomNumGen->randNorm(0, random_component) ,
				randomNumGen->randNorm(0, random_component) ,
				randomNumGen->randNorm(0, random_component));
		protein->set_bb_atom_coords(pos+randVec,POSE::CA, r);
		pos = pos + incr;
	}

	vector3d com = get_ca_CoM(protein);
	for ( int r = 0; r < protein->get_residue_count(); r++){
		protein->get_bb_atom(POSE::CA, r)->set_coords(protein->get_bb_atom(POSE::CA, r)->get_coords() - com );
	}
	bool_vector atom_selection(protein->get_all_atom_count(), true);
	LOOP::ca_minimise_bonds_bumps_angles(protein, atom_selection);
}

std::ostream& output_psicov_contact_map(PRODART::POSE::const_pose_shared_ptr protein,
		const double cb_cutoff,
		const int min_seq_sep,
		std::ostream& output){

	protein->index();

	for (int i = 0; i < protein->get_residue_count(); i++){
		const_residue_shared_ptr res_i = protein->get_residue(i);
		const_atom_shared_ptr cb_i = res_i->get_type().is_equal3(residue_type("GLY")) ? res_i->get_bb_atom(POSE::CA) : res_i->get_bb_atom(POSE::CB);

		if (cb_i->isActiveAndSet()){
			for (int j = i+min_seq_sep; j < protein->get_residue_count(); j++){
				const_residue_shared_ptr res_j = protein->get_residue(j);
				const_atom_shared_ptr cb_j = res_j->get_type().is_equal3(residue_type("GLY")) ? res_j->get_bb_atom(POSE::CA) : res_j->get_bb_atom(POSE::CB);
				if (cb_j->isActiveAndSet()){
					const double dist = (cb_i->get_coords() - cb_j->get_coords()).mod();
					if (dist < cb_cutoff){
						output  << res_i->get_internal_residue_index() + 1 << " "
								<< res_j->get_internal_residue_index() + 1 << " "
								<< "0 "
								<< cb_cutoff << " "
								<< "1.0\n";
					}
				}
			}
		}
	}

	return output;
}


std::map< int,std::map< int, boost::tuple<double, double, double> > > make_psicov_map(PRODART::POSE::const_pose_shared_ptr protein,
		const double cb_cutoff,
		const int min_seq_sep){
	std::map< int,std::map< int, boost::tuple<double, double, double> > > rtn_map;

	protein->index();

	for (int i = 0; i < protein->get_residue_count(); i++){
		const_residue_shared_ptr res_i = protein->get_residue(i);
		const_atom_shared_ptr cb_i = res_i->get_type().is_equal3(residue_type("GLY")) ? res_i->get_bb_atom(POSE::CA) : res_i->get_bb_atom(POSE::CB);

		if (cb_i->isActiveAndSet()){
			for (int j = i+min_seq_sep; j < protein->get_residue_count(); j++){
				const_residue_shared_ptr res_j = protein->get_residue(j);
				const_atom_shared_ptr cb_j = res_j->get_type().is_equal3(residue_type("GLY")) ? res_j->get_bb_atom(POSE::CA) : res_j->get_bb_atom(POSE::CB);
				if (cb_j->isActiveAndSet()){
					const double dist = (cb_i->get_coords() - cb_j->get_coords()).mod();
					if (dist < cb_cutoff){
						const int resno_i = lexical_cast<int>(res_i->get_trimmed_pdb_residue_index());
						const int resno_j = lexical_cast<int>(res_j->get_trimmed_pdb_residue_index());
						const double lower = 0;
						const double upper = cb_cutoff;
						const double prob = 1.0;
		                boost::tuple<double, double, double> tup =  boost::make_tuple(lower, upper, prob);
		                rtn_map[resno_i][resno_j] = tup;
					}
				}
			}
		}
	}



	return rtn_map;
}

std::map< int,std::map< int, boost::tuple<double, double, double> > > parse_psicov_file(std::istream& psicov_input){
	std::map< int,std::map< int, boost::tuple<double, double, double> > > rtn_map;

    string lineStr;
    unsigned long length, lineNum = 0 ;
    string_vector SplitVec;
    //int pdbs_loaded = 0;

	//cout << "reading file..." << endl;
    while ( !psicov_input.eof() ) {
            getline(psicov_input, lineStr);

            trim(lineStr);
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {

                split( SplitVec, lineStr, is_any_of("\t ") );

                int res1 = lexical_cast<int>(SplitVec[0]);
                int res2 = lexical_cast<int>(SplitVec[1]);
                double lower = lexical_cast<double>(SplitVec[2]);
                double upper = lexical_cast<double>(SplitVec[3]);
                double prob = lexical_cast<double>(SplitVec[4]);

                boost::tuple<double, double, double> tup =  boost::make_tuple(lower, upper, prob);

                rtn_map[res1][res2] = tup;

            }
    }
	return rtn_map;
}

std::ostream& compare_to_psicov_contact_map(PRODART::POSE::const_pose_shared_ptr protein,
		const double cb_cutoff,
		const double prob_cutoff,
		const int min_seq_sep,
		std::istream& psicov_input,
		std::ostream& output){

	std::map< int,std::map< int, boost::tuple<double, double, double> > > psicov_mat = parse_psicov_file(psicov_input);
	std::map< int,std::map< int, boost::tuple<double, double, double> > > real_mat = make_psicov_map(protein, cb_cutoff, min_seq_sep);

	typedef std::map< int,std::map< int, boost::tuple<double, double, double> > >::iterator outer_it;
	typedef std::map< int, boost::tuple<double, double, double> >::iterator inner_it;

	unsigned int false_pos_count = 0, false_neg_count = 0, true_pos_count = 0, true_neg_count = 0;

	unsigned int total_psicov = 0, total_real = 0;


	for (outer_it it_i = psicov_mat.begin(); it_i != psicov_mat.end(); it_i++){
		const int res_i = it_i->first;
		for (inner_it it_j = it_i->second.begin(); it_j != it_i->second.end(); it_j++){
			const int res_j = it_j->first;
			const int seq_sep = res_j - res_i;

			if (seq_sep >= min_seq_sep){
				const double prob = it_j->second.get<2>();
				if (prob >= prob_cutoff){
					total_psicov++;
				}
				//does real_mat have this contact?
				if (real_mat.find(res_i) != real_mat.end()){
					if (real_mat.find(res_i)->second.find(res_j) != real_mat.find(res_i)->second.end()){
						if (prob >= prob_cutoff){
							// in both mats
							true_pos_count++;
						}
						else {
							if (prob >= prob_cutoff) false_pos_count++;
						}
					}
					else {
						if (prob >= prob_cutoff) false_pos_count++;
					}
				}
				else {
					if (prob >= prob_cutoff) false_pos_count++;
				}
			}
		}
	}


	for (int it_i = 0; it_i < protein->get_residue_count(); it_i++){
		const int res_i = lexical_cast<int>(protein->get_residue(it_i)->get_trimmed_pdb_residue_index());
		for (int it_j = it_i + 1; it_j != protein->get_residue_count(); it_j++){
			const int res_j = lexical_cast<int>(protein->get_residue(it_j)->get_trimmed_pdb_residue_index());
			const int seq_sep = res_j - res_i;
			if (seq_sep >= min_seq_sep){
				bool is_contact = false;
				if (real_mat.find(res_i) != real_mat.end()){
					if (real_mat.find(res_i)->second.find(res_j) !=  real_mat.find(res_i)->second.end()){
						is_contact = true;
						total_real++;
					}
				}

				//does psicov_mat have this contact?
				if (psicov_mat.find(res_i) != psicov_mat.end()){
					if (psicov_mat.find(res_i)->second.find(res_j) != psicov_mat.find(res_i)->second.end()){
						const double prob = psicov_mat.find(res_i)->second.find(res_j)->second.get<2>();
						if (prob >= prob_cutoff){
							if (is_contact ==  false){
								//false in real; true in psicov
								//false_pos_count++;
							}
							else {
								// true in real, true in psicov
								//true_pos_count++;
							}
						}
						else {
							if (is_contact ==  false){
								//false in real; false in psicov
								true_neg_count++;
							}
							else {
								// true in real, false in psicov
								false_neg_count++;
							}
						}
					}
					else {
						if (is_contact ==  false){
							//false in real; false in psicov
							true_neg_count++;
						}
						else {
							// true in real, false in psicov
							false_neg_count++;
						}
					}
				}
				else {
					if (is_contact ==  false){
						//false in real; false in psicov
						true_neg_count++;
					}
					else {
						// true in real, false in psicov
						false_neg_count++;
					}
				}
			}
		}
	}


	output << "prob_cutoff\ttrue_pos\ttrue_neg\tfalse_pos\tfalse_neg\ttotal_psicov\ttotal_real\n";
	output << prob_cutoff << "\t"
			<< double(true_pos_count) /  double(1) << "\t"
			<< double(true_neg_count) /  double(1) << "\t"
			<< double(false_pos_count) / double(1) << "\t"
			<< double(false_neg_count) / double(1) << "\t"
			<< total_psicov << "\t"
			<< total_real << "\t"
			<< endl;

	return output;
}

std::ostream& output_contactmap_rstfile(PRODART::POSE::pose_shared_ptr protein,
		const double cutoff,
		std::ostream& output){
	protein->index();

	for (int res_i= 0; res_i < protein->get_residue_count(); res_i++){
		atom_shared_ptr at_i = protein->get_bb_atom(POSE::CA, res_i);
		for (int res_j = res_i+1; res_j < protein->get_residue_count(); res_j++){
			atom_shared_ptr at_j = protein->get_bb_atom(POSE::CA, res_j);
			if ((res_j - res_i) > 4){
				const double dist = (at_i->get_coords() - at_j->get_coords()).mod();
				if (dist < cutoff){
					output << "CA_GO_CONTACT "
							<< at_i->get_residue()->get_trimmed_pdb_residue_index() << " " << at_i->get_chain()->getChainID() << " "
							<< at_j->get_residue()->get_trimmed_pdb_residue_index() << " " << at_j->get_chain()->getChainID() << " "
							<< "0.0 "
							<< cutoff << " "
							<< "1.0 "
							<< endl;
				}
			}
		}
	}

	return output;
}


PRODART::POSE::pose_shared_ptr make_dummy_ca_pose(const PRODART::POSE::residue_type_vector seq,
		const char chainID,
		const int start_res,
		const double rand_coord){


    pose_shared_ptr pose_ = PRODART::POSE::new_pose();
    pose_->add_new_chain(chainID);


    for (unsigned int i = 0; i < seq.size(); i++){
            residue_shared_ptr res = pose_->append_residue(seq[i], chainID);
            pose_->index();
            vector3d rand(ENV::get_random_num_gen()->rand(rand_coord),ENV::get_random_num_gen()->rand(rand_coord),ENV::get_random_num_gen()->rand(rand_coord));
            pose_->add_new_atom(rand, atom_type("N"), res->get_internal_residue_index());
            rand = vector3d(ENV::get_random_num_gen()->rand(rand_coord),ENV::get_random_num_gen()->rand(rand_coord),ENV::get_random_num_gen()->rand(rand_coord));
            pose_->add_new_atom(rand, atom_type("CA"), res->get_internal_residue_index());
            rand = vector3d(ENV::get_random_num_gen()->rand(rand_coord),ENV::get_random_num_gen()->rand(rand_coord),ENV::get_random_num_gen()->rand(rand_coord));
            pose_->add_new_atom(rand, atom_type("C"), res->get_internal_residue_index());
            pose_->index();
    }

    pose_->renumber_residues(chainID, start_res);
    pose_->index();

    return pose_;
}


bool thread_sequence(const PRODART::POSE::residue_type_vector seq,
		PRODART::POSE::pose_shared_ptr pose_){

	if (seq.size() != pose_->get_residue_count()){
		cerr << "thread_sequence: ERROR: sequence not the same size as the pose. seq_len:" << seq.size()
				<< " pose_len:" << pose_->get_residue_count()
				<< endl;
		return false;
	}

	for (int i = 0; i < pose_->get_residue_count(); i++){
		if (!pose_->get_residue(i)->get_type().is_equal3(seq[i])){
			pose_->get_residue(i)->get_sidechain()->clear();
			pose_->get_residue(i)->set_type(seq[i]);
			if (pose_->get_residue(i)->get_type().is_equal3(residue_type("GLY"))){
				pose_->get_bb_atom(POSE::CB, i)->setActive(false);
				pose_->get_bb_atom(POSE::CB, i)->setSet(false);
			}
		}
	}

	return true;
}


double get_rmsd_superpose_atmmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::pose_shared_ptr protein_to_move,
		std::istream& input_atom_mapping){
	std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;

    string lineStr;

    //long length = 0;//, lineNum = 0 ;


	while ( !input_atom_mapping.eof() ) {
		getline(input_atom_mapping, lineStr);

		//length = lineStr.length();

		string_vector SplitVec;
		split( SplitVec, lineStr, is_any_of("\t ") );

		if (SplitVec.size() >= 6){
			if ( SplitVec[0].substr(0,1).compare("#") != 0){
				const string ref_resnum = SplitVec[0];
				char ref_chain = lexical_cast<char>(SplitVec[1]);
				const atom_type ref_atmStr(SplitVec[2]);
				if (ref_chain == '-') ref_chain = ' ';
				const string resnum = SplitVec[3];
				char tm_chain = lexical_cast<char>(SplitVec[4]);
				if (tm_chain == '-') tm_chain = ' ';
				const atom_type tm_atmStr(SplitVec[5]);

				const_atom_shared_ptr ref_atom = ref_protein->get_atom(ref_atmStr, ref_resnum, ref_chain);
				const_atom_shared_ptr tm_atom = protein_to_move->get_atom(tm_atmStr, resnum, tm_chain);

				if (ref_atom && tm_atom){
					if (ref_atom->isSet() && tm_atom->isSet()){
						atom_mapping[ref_atom] = tm_atom;
						cout << "mapped:\t"
								<< ref_atom->get_residue()->get_type().get_label() << "\t"
								<< ref_atom->get_residue()->get_pdb_residue_index() << "\t"
								<< ref_chain << "\t"
								<< ref_atom->get_type().get_pdb_formatted_label() << "\t"
								<< "->" << "\t"
								<< tm_atom->get_residue()->get_type().get_label() << "\t"
								<< tm_atom->get_residue()->get_pdb_residue_index() << "\t"
								<< tm_chain << "\t"
								<< tm_atom->get_type().get_pdb_formatted_label()
								<< endl;
					}

				}
				else {
					cerr << "ERROR: get_rmsd_superpose: residues not found\n";

					if (!ref_atom){
						cerr << "\t"
							 << "'" << ref_resnum << "'\t"
							 << "'" << ref_chain << "'\t"
							 << "'" << ref_atmStr.get_pdb_formatted_label() << "'\t"
							 << "in reference pose: " << ref_protein->get_label() << "\t"
							 << "not found\n";
					}

					if (!tm_atom){
						cerr << "\t"
							 << "'" << resnum << "'\t"
							 << "'" << tm_chain << "'\t"
							 << "'" << tm_atmStr.get_pdb_formatted_label() << "'\t"
							 << "in move pose: " << protein_to_move->get_label() << "\t"
							 << "not found\n";
					}

				}

			}

		}

	}

	return get_rmsd_superpose(protein_to_move, atom_mapping);
}

double get_rmsd_superpose_resmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::pose_shared_ptr protein_to_move,
		std::istream& input_residue_mapping,
		bool verbose){


    std::vector<boost::tuple< std::string, char, std::string, char  > > res_map = parse_residue_mapping(input_residue_mapping);
    std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping = get_const_atom_mapping(ref_protein, protein_to_move, res_map, verbose);

	return get_rmsd_superpose(protein_to_move, atom_mapping);

}

double get_rmsd_no_superpose_resmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::const_pose_shared_ptr protein_to_move,
		std::istream& input_residue_mapping,
		bool verbose){

    std::vector<boost::tuple< std::string, char, std::string, char  > > res_map = parse_residue_mapping(input_residue_mapping);
    std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping = get_const_atom_mapping(ref_protein, protein_to_move, res_map, verbose);

	return get_fixed_position_rmsd(atom_mapping);
}

std::vector<boost::tuple< std::string, char, std::string, char  > > parse_residue_mapping(std::istream& input_residue_mapping){
	std::vector<boost::tuple< std::string, char, std::string, char  > > rtn_map;

    string lineStr;
	while ( !input_residue_mapping.eof() ) {
		getline(input_residue_mapping, lineStr);

		//length = lineStr.length();

		string_vector SplitVec;
		split( SplitVec, lineStr, is_any_of("\t ") );

		if (SplitVec.size() >= 4){
			if ( SplitVec[0].substr(0,1).compare("#") != 0){
				const string ref_resnum = SplitVec[0];
				char ref_chain = lexical_cast<char>(SplitVec[1]);
				if (ref_chain == '-') ref_chain = ' ';
				const string resnum = SplitVec[2];
				char tm_chain = lexical_cast<char>(SplitVec[3]);
				if (tm_chain == '-') tm_chain = ' ';

				boost::tuple< std::string, char, std::string, char  > tup = boost::tuple< std::string, char, std::string, char  > (ref_resnum, ref_chain, resnum, tm_chain  );

				rtn_map.push_back(tup);
			}
		}
	}

	return rtn_map;
}

std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> get_const_atom_mapping(
		const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::const_pose_shared_ptr protein_to_move,
		std::vector<boost::tuple< std::string, char, std::string, char  > > residue_map,
		bool verbose){
	std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;

	std::vector<boost::tuple< std::string, char, std::string, char  > >::const_iterator iter;

	for (iter = residue_map.begin(); iter != residue_map.end(); iter++){
		const string ref_resnum = iter->get<0>();
		char ref_chain = iter->get<1>();
		if (ref_chain == '-') ref_chain = ' ';
		const string resnum = iter->get<2>();
		char tm_chain = iter->get<3>();
		if (tm_chain == '-') tm_chain = ' ';

		if (ref_protein->get_residue(ref_resnum, ref_chain)
				&& protein_to_move->get_residue(resnum, tm_chain)){
			const int start_bb_at = ref_protein->get_residue(ref_resnum, ref_chain)->get_first_bb_atom_index();
			for (int bb_atom_num = start_bb_at; bb_atom_num < start_bb_at + pose::get_num_protein_backbone_atom_types(); bb_atom_num++ ){
				const_atom_shared_ptr ref_atom = ref_protein->get_bb_atom(bb_atom_num);
				const_atom_shared_ptr tm_atom = protein_to_move->get_atom(ref_atom->get_type(),resnum, tm_chain);
				if (ref_atom && tm_atom){
					if (ref_atom->isSet() && tm_atom->isSet()){
						atom_mapping[ref_atom] = tm_atom;
						if (verbose){
							cout << "mapped:\t"
									<< ref_atom->get_residue()->get_type().get_label() << "\t"
									<< ref_atom->get_residue()->get_pdb_residue_index() << "\t"
									<< ref_chain << "\t"
									<< ref_atom->get_type().get_pdb_formatted_label() << "\t"
									<< "->" << "\t"
									<< tm_atom->get_residue()->get_type().get_label() << "\t"
									<< tm_atom->get_residue()->get_pdb_residue_index() << "\t"
									<< tm_chain << "\t"
									<< tm_atom->get_type().get_pdb_formatted_label()
									<< endl;
						}
					}

				}
			}
			//  add in sidechain atom mapping here
			const int num_sc_atoms = ref_protein->get_residue(ref_resnum, ref_chain)->get_sidechain()->get_atom_count();
			for (int i = 0; i < num_sc_atoms; i++){
				const_atom_shared_ptr ref_atom = ref_protein->get_residue(ref_resnum, ref_chain)->get_sidechain()->get_atom(i);
				const_atom_shared_ptr tm_atom = protein_to_move->get_atom(ref_atom->get_type(),resnum, tm_chain);
				if (ref_atom && tm_atom){
					if (ref_atom->isSet() && tm_atom->isSet()){
						atom_mapping[ref_atom] = tm_atom;
						if (verbose){
							cout << "mapped:\t"
									<< ref_atom->get_residue()->get_type().get_label() << "\t"
									<< ref_atom->get_residue()->get_pdb_residue_index() << "\t"
									<< ref_chain << "\t"
									<< ref_atom->get_type().get_pdb_formatted_label() << "\t"
									<< "->" << "\t"
									<< tm_atom->get_residue()->get_type().get_label() << "\t"
									<< tm_atom->get_residue()->get_pdb_residue_index() << "\t"
									<< tm_chain << "\t"
									<< tm_atom->get_type().get_pdb_formatted_label()
									<< endl;
						}
					}

				}

			}

		}
		else {
			cerr << "ERROR: get_rmsd_superpose: residues not found\n";

			if (!ref_protein->get_residue(ref_resnum, ref_chain)){
				cerr << "\t"
						<< "'" << ref_resnum << "'\t"
						<< "'" << ref_chain << "'\t"
						<< "in reference pose: " << ref_protein->get_label() << "\t"
						<< "not found\n";
			}

			if (!protein_to_move->get_residue(resnum, tm_chain)){
				cerr << "\t"
						<< "'" << resnum << "'\t"
						<< "'" << tm_chain << "'\t"
						<< "in move pose: " << protein_to_move->get_label() << "\t"
						<< "not found\n";
			}

		}
		// end else block

	}

	return atom_mapping;
}


std::vector<double> get_rmsd_no_superpose_trajectory_resmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		std::istream& traj_input,
		std::istream& input_residue_mapping,
		bool verbose){
	std::vector<double> rmsds_vec;

    std::vector<boost::tuple< std::string, char, std::string, char  > > res_map = parse_residue_mapping(input_residue_mapping);

    unsigned long count = 1;

	while ( !traj_input.eof()){



		pose_shared_ptr pose_ = new_pose();
		pose_->loadPdb(traj_input);

		if (pose_->get_residue_count() > 0){
			std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping = get_const_atom_mapping(ref_protein, pose_, res_map, verbose);
			double rmsd = get_fixed_position_rmsd(atom_mapping);

			rmsds_vec.push_back(rmsd);

			cout << "traj_RMSD:\t" << count << "\t" << rmsd << endl;
		}

	    count++;

	}

	return rmsds_vec;
}


std::ostream& output_rosetta_distance_restraints_all_atom(PRODART::POSE::const_pose_shared_ptr protein,
		const double max_dist_cutoff,
		const double stddev,
		std::ostream& output){


	PRODART::POSE::pose_shared_ptr prot_clone = protein->clone();
	POSE_UTILS::quick_add_HN(prot_clone,false);
	for (int i = 0; i < prot_clone->get_residue_count(); i++){
		atom_shared_ptr at = prot_clone->get_bb_atom(POSE::H, i);
		at->set_type(atom_type("H"));
	}

	bb_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_bb_pose_meta(prot_clone);
	PRODART::POSE::three_state_sec_struct_vector secs = pose_meta_->get_sec_struct();

	const unsigned long atom_count = prot_clone->get_all_atom_count();

	for (unsigned long ii = 0; ii < atom_count; ii++){
		const_atom_shared_ptr atom_i = prot_clone->get_atom(ii);
		if (atom_i->isActiveAndSet()){
			const int res_i = atom_i->get_residue()->get_internal_residue_index();
			const vector3d pos_i = atom_i->get_coords();
			for (unsigned long jj = ii+1; jj < atom_count; jj++){
				const_atom_shared_ptr atom_j = prot_clone->get_atom(jj);
				const int res_j = atom_j->get_residue()->get_internal_residue_index();
				if (atom_j->isActiveAndSet() && std::abs(res_j - res_i) >= 4){
					const vector3d pos_j = atom_j->get_coords();
					const double dist = (pos_i - pos_j).mod();
					if (dist <= max_dist_cutoff){

						if
							( ( atom_i->get_type().get_trimmed_label().find("H") == string::npos
								&& atom_j->get_type().get_trimmed_label().find("H") == string::npos )
								|| (atom_i->get_type().get_trimmed_label().compare("H") == 0 && atom_j->get_type().get_trimmed_label().find("H") == string::npos )
								|| (atom_j->get_type().get_trimmed_label().compare("H") == 0 && atom_i->get_type().get_trimmed_label().find("H") == string::npos ) )

							{
							output << "AtomPair\t"
									<< atom_i->get_type().get_trimmed_label() << "\t" << atom_i->get_residue()->get_pdb_residue_index() << atom_i->get_chain()->getChainID() << "\t"
									<< atom_j->get_type().get_trimmed_label() << "\t" << atom_j->get_residue()->get_pdb_residue_index() << atom_j->get_chain()->getChainID() << "\t"
									<< "HARMONIC\t" << dist << "\t" << stddev
									<< "\t#"
									<< atom_i->get_residue()->get_type().get_label3()
									//<< atom_i->get_chain()->getChainID()
									<< "\t"
									<< atom_j->get_residue()->get_type().get_label3()
									//<< atom_j->get_chain()->getChainID()
									<< "\t";

							if (secs[res_j] == secs[res_i] ){

								if (secs[res_i] == ss3_HELIX ){
									output << "# HELIX";
								}
								else if (secs[res_i] == ss3_STRAND) {
									output << "# STRAND";
								}
								else  {
									output << "# OTHER";
								}
							}
							else {
								output << "# OTHER";
							}

							output << std::endl;
						}
					}
				}
			}
		}
	}

	return output;

}


std::ostream& output_rosetta_coord_restraints_all_atom(PRODART::POSE::const_pose_shared_ptr protein,
		const double stddev,
		std::ostream& output){
	const unsigned long atom_count = protein->get_all_atom_count();

	const_residue_shared_ptr first_residue = protein->get_residue(0);
	const_residue_shared_ptr second_residue = protein->get_residue(1);
	//protein->get_residue(0)->get_bb_atom(BBAtomType::CA);


	for (unsigned long ii = 0; ii < atom_count; ii++){
		const_atom_shared_ptr atom_i = protein->get_atom(ii);
		if (atom_i->isActiveAndSet()){
			const vector3d pos_i = atom_i->get_coords();

			if (atom_i->get_type().get_trimmed_label().find("H") == string::npos){

				string ref_resid = atom_i->get_residue() != first_residue ? (first_residue->get_pdb_residue_index() + first_residue->get_chain()->getChainID()) : (second_residue->get_pdb_residue_index() + second_residue->get_chain()->getChainID() )  ;

				output << "CoordinateConstraint\t"
						<< atom_i->get_type().get_trimmed_label() << "\t" << atom_i->get_residue()->get_pdb_residue_index()
						<< atom_i->get_chain()->getChainID()
						<< "\t"
						<< "CA\t" << ref_resid << "\t"
						<< pos_i.x << "\t" << pos_i.y << "\t" << pos_i.z << "\t"
						<< "HARMONIC\t" << 0 << "\t" << stddev
						<< "\t#"
						<< atom_i->get_residue()->get_type().get_label3() << "\t"
						<< endl;
			}

		}
	}

	return output;

}






bool search_replace_atom_types(const PRODART::POSE::pose_shared_ptr pose_, std::istream& search_replace_mapping){

    string lineStr;


    unsigned long length, lineNum = 0 ;

    std::map<POSE::atom_type, POSE::atom_type> s_r_map;


    string_vector SplitVec;
	cout << "reading file..." << endl;
    while ( !search_replace_mapping.eof() ) {
            getline(search_replace_mapping, lineStr);
            string resStr;
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {
                    split( SplitVec, lineStr, is_any_of("\t") );
                    if ( SplitVec[0].substr(0,1).compare("#") != 0
                     && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 2 ){

                    	atom_type at_type = atom_type(SplitVec[0]);
                    	atom_type rep_at_type = atom_type(SplitVec[1]);

                    	if (s_r_map.find(at_type) != s_r_map.end()){
                    		cerr << "misc_protocols::search_replace_atom_types: WARNING: you have duplicate search terms. This is a duplicate search atom_type: "
                    				<< at_type.get_trimmed_label() << "->" << rep_at_type.get_trimmed_label()
                    				<< endl;
                    	}

                    	s_r_map[at_type] = rep_at_type;

                    }
            }
    }

    cout << "misc_protocols::search_replace_atom_types: carrying out these search/replace operations: " << endl;
    std::map<POSE::atom_type, POSE::atom_type>::const_iterator iter;
    for (iter = s_r_map.begin(); iter != s_r_map.end(); iter++){
    	cout << iter->first.get_trimmed_label() << "->" << iter->second.get_trimmed_label() << endl;
    }

    bool result = POSE_UTILS::search_replace_atom_types(pose_, s_r_map);

    if (result == false){
    	cerr << "misc_protocols::search_replace_atom_types: WARNING: there were problems with this search/replace operation." << endl;
    }

    return result;

}












}
}
}

