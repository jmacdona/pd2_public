/*
 * alphabet_bb_builder.cpp
 *
 *  Created on: May 2, 2012
 *      Author: jmacdona
 */
#include "alphabet_bb_builder.h"



using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;

namespace PRODART {
namespace POSE {
namespace BB_BUILDER {



boost::shared_ptr<alphabet_bb_builder> new_alphabet_bb_builder(){
	boost::shared_ptr<alphabet_bb_builder> nalphabet_bb_builder(new alphabet_bb_builder);
	//nalphabet_bb_builder->_this = nalphabet_bb_builder;
	return nalphabet_bb_builder;
}


alphabet_bb_builder::alphabet_bb_builder(){
	pose_alphabet.clear();
	letter_length = 0;
	centre_weight = 10;
	is_even = false;
	default_b_factor = 0.0;
}

class load_alphabet_exception : public std::exception{
  virtual const char* what() const throw()
  {
    return "ERROR: PDB alphabet file contained errors";
  }
};

bool alphabet_bb_builder::load_alphabet_pose(std::istream & input){
	pose_alphabet.clear();

	int letter_count = 0;
	int frag_size = 0;

	while ( !input.eof() ) {
		pose_shared_ptr letter = new_pose();
		letter->set_label("TEMP_NAME");
		letter->loadPdb(input);
		//letter->set_label("TEMP_NAME2");
		if (letter->get_residue_count() != 0){
			//cout << "DEBUG1" << endl;
			pose_alphabet.push_back(letter);
			//cout << "DEBUG2" << endl;
			//PRINT_EXPR(letter->get_bb_atom(POSE::CA, 3)->get_coords());
			letter_count++;
			if (frag_size != letter->get_residue_count() && frag_size != 0){
				cerr << "\nalphabet_bb_builder: WARNING: alphabet seems to contain letters of different sizes!!!\n" << endl;
				throw load_alphabet_exception();
			}
			frag_size = letter->get_residue_count();
			std::stringstream lab;
			lab << "Letter:" << letter_count << " Letter_size:" << frag_size;
			letter->set_label(lab.str());
		}
	}

	letter_length = frag_size;

	is_even = letter_length % 2 == 0 ? true : false;
	const int centre_res = letter_length / 2;

	weights = double_vector(letter_length, 1);

	if (is_even){
		weights[centre_res-1] = centre_weight;
		weights[centre_res] = centre_weight;
	}
	else {
		weights[centre_res-1] = centre_weight;
		weights[centre_res] = centre_weight;
		weights[centre_res+1] = centre_weight;
	}


	double curr_weight = 1;
	for (int i = 2; i < letter_length; i++){
		curr_weight = curr_weight / centre_weight;
		if (is_even){
			const int lower = centre_res - i - 1;
			const int upper = centre_res + i;

			if (lower < 0 || upper >= letter_length ){
				break;
			}

			weights[lower] = curr_weight;
			weights[upper] = curr_weight;

		}
		else {
			const int lower = centre_res - i;
			const int upper = centre_res + i;
			if (lower < 0 || upper >= letter_length ){
				break;
			}

			weights[lower] = curr_weight;
			weights[upper] = curr_weight;

		}

	}

	/*
	for (int i = 0; i < letter_length; i++){
		PRINT_EXPR(weights[i]);
	}
	*/


	cout << "alphabet_bb_builder: loaded alphabet of " << letter_count << " letters of length " << frag_size << endl;

	return true;
}

double alphabet_bb_builder::get_rmsd(const int first_res, const_pose_shared_ptr letter, pose_shared_ptr pose_) const{
	if (first_res + letter_length > pose_->get_residue_count()){
		cerr << "alphabet_bb_builder::get_rmsd: ERROR: out of range" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if (pose_->get_residue(first_res)->get_chain() != pose_->get_residue(first_res + letter_length -1)->get_chain()){
		cerr << "alphabet_bb_builder::get_rmsd: ERROR: residues on different chains" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	//const bool is_even = letter_length % 2 == 0 ? true : false;
	const int centre_res = letter_length / 2;

	vector<boost::tuple<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr, double> > at_at_wt_vec;
	at_at_wt_vec.reserve(letter_length);
	//std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;
	for (int i = 0; i < letter_length; i++){
		const int resnum = first_res + i;
		boost::tuple<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr, double> tup(letter->get_bb_atom(POSE::CA, i),
				pose_->get_bb_atom(POSE::CA, resnum),
				weights[i]);
		at_at_wt_vec.push_back(tup);
		//atom_mapping[letter->get_bb_atom(POSE::CA, i)] = pose_->get_bb_atom(POSE::CA, resnum);
	}



	const vector3d CoM1 = is_even ?  letter->get_bb_atom(POSE::CA, centre_res-1)->get_coords() :  letter->get_bb_atom(POSE::CA, centre_res)->get_coords();
	const vector3d CoM2 = is_even ?  pose_->get_bb_atom(POSE::CA, first_res + centre_res-1)->get_coords() :  pose_->get_bb_atom(POSE::CA, first_res + centre_res)->get_coords();


	return sqrt(POSE_UTILS::get_msd(at_at_wt_vec, CoM1, CoM2));

}

pose_shared_ptr alphabet_bb_builder::get_best_fit_letter(const int first_res, pose_shared_ptr pose_) const{
	double best_rmsd = std::numeric_limits<double>::max();
	unsigned int best_letter = 0;
	for (unsigned int i = 0; i < pose_alphabet.size(); i++){
		const double rmsd = this->get_rmsd(first_res, pose_alphabet[i], pose_);
		/*
		cout << "DEBUGGER\t" << first_res << "\t" << i << "\t" << rmsd << "\t"
				<< pose_alphabet[i]->get_bb_atom(POSE::CA,0)->get_coords() << "\t"
				<< pose_->get_bb_atom(POSE::CA,first_res)->get_coords() << "\t"
				<< endl;
				*/
		if (rmsd < best_rmsd){
			best_rmsd = rmsd;
			best_letter = i;
		}
	}
	//PRINT_EXPR(best_letter);
	//PRINT_EXPR(best_rmsd);
	//PRINT_EXPR(pose_alphabet[best_letter]->get_bb_atom(POSE::CA, 3)->get_coords());
	//pose_alphabet[best_letter]->outputPdb(cout);
	pose_shared_ptr b_let = pose_alphabet[best_letter]->clone();
	//PRINT_EXPR(b_let->get_bb_atom(POSE::CA, 3)->get_coords());

	//PRINT_EXPR(b_let->get_residue_count());

	//const bool is_even = letter_length % 2 == 0 ? true : false;
	const int centre_res = letter_length / 2;

	vector<boost::tuple<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr, double> > at_at_wt_vec;
	at_at_wt_vec.reserve(letter_length);

	//std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;
	for (int i = 0; i < letter_length; i++){
		const int resnum = first_res + i;
		boost::tuple<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr, double> tup(pose_->get_bb_atom(POSE::CA, resnum),
				b_let->get_bb_atom(POSE::CA, i),
				weights[i]);
		at_at_wt_vec.push_back(tup);
		//atom_mapping[pose_->get_bb_atom(POSE::CA, resnum)] = b_let->get_bb_atom(POSE::CA, i);
	}



	const vector3d CoM_let = is_even ?  b_let->get_bb_atom(POSE::CA, centre_res-1)->get_coords() :  b_let->get_bb_atom(POSE::CA, centre_res)->get_coords();
	const vector3d CoM_pose = is_even ?  pose_->get_bb_atom(POSE::CA, first_res + centre_res-1)->get_coords() :  pose_->get_bb_atom(POSE::CA, first_res + centre_res)->get_coords();

	//b_let->outputPdb(cout);
	//const double rmsd =
	get_rmsd_superpose(b_let, at_at_wt_vec, CoM_pose, CoM_let);
	//PRINT_EXPR(rmsd);
	//PRINT_EXPR(b_let->get_residue_count());

	return b_let;

}

POSE::atom_shared_ptr_vector alphabet_bb_builder::build_fragment(const int first_res, pose_shared_ptr pose_) const{

	atom_shared_ptr_vector rtn_vec;
	rtn_vec.reserve(10);

	//PRINT_EXPR("");

	pose_shared_ptr letter = this->get_best_fit_letter(first_res, pose_);

	//PRINT_EXPR(letter->get_residue_count());

	const int last_res =  first_res + letter_length - 1;

	if (pose_->get_residue(first_res)->get_chain() != pose_->get_residue(last_res)->get_chain()){
		//cerr << "alphabet_bb_builder::build_fragment: ERROR: attemping to build across different chains" << endl;
		return rtn_vec;
	}

	for (int i =0; i < letter_length; i++){
		const int pose_res_num = i + first_res;
		//const int last_res = i + first_ca + letter_length - 1;
		residue_shared_ptr res = letter->get_residue(i);
		residue_shared_ptr pose_res = pose_->get_residue(pose_res_num);

		atom_shared_ptr_vector atoms_to_add = res->get_all_atoms();
		for (atom_shared_ptr_vector::iterator iter = atoms_to_add.begin(); iter != atoms_to_add.end(); iter++){
			atom_shared_ptr atm = *iter;
			if (atm->isActiveAndSet()){
				if (atm->get_type() != atom_type("CA")){
					if (!((pose_res->get_type() == residue_type("GLY")) && atm->get_type() == atom_type("CB"))){
						atom_shared_ptr nat = pose_->add_new_atom(atm->get_coords(), atm->get_type(), pose_res_num);
						nat->set_b_factor(default_b_factor);
						rtn_vec.push_back(nat);
					}
				}

			}
		}



	}
	return rtn_vec;
}

//! returns rebuild atoms vector
POSE::atom_shared_ptr_vector alphabet_bb_builder::build_peptide_bond(const int first_res, pose_shared_ptr pose_) const{
	if (this->letter_length % 2 != 0){
		cerr << "alphabet_bb_builder::build_peptide_bond: ERROR: need an even numbered alphabet" << endl;
		return atom_shared_ptr_vector();
	}
	const int real_first_res = first_res - ((letter_length/2)-1);
	return this->build_fragment(real_first_res, pose_);
}

POSE::atom_shared_ptr_vector alphabet_bb_builder::build_peptide_bonds(const int first_res, const int last_res, pose_shared_ptr pose_) const{
	atom_shared_ptr_vector rtn_vec;
	rtn_vec.reserve(10 * (last_res-first_res));

	if (this->letter_length % 2 != 0){
		cerr << "alphabet_bb_builder::build_peptide_bonds: ERROR: need to use an even numbered alphabet" << endl;
		return rtn_vec;
	}
	else if (last_res <= first_res){
		cerr << "alphabet_bb_builder::build_peptide_bonds: ERROR: last_res must be more than first_res" << endl;
		return rtn_vec;
	}
	else if (last_res >= pose_->get_residue_count()){
		cerr << "alphabet_bb_builder::build_peptide_bonds: ERROR: last_res out of range" << endl;
		return rtn_vec;
	}
	else if (first_res < 0){
		cerr << "alphabet_bb_builder::build_peptide_bonds: ERROR: first_res out of range" << endl;
		return rtn_vec;
	}



	if (pose_->get_residue(first_res)->get_chain() != pose_->get_residue(last_res)->get_chain()){
		cerr << "alphabet_bb_builder::build_peptide_bonds: WARNING: first_res and last_res are on different chains" << endl;
	}

	for (int i = first_res; i < last_res; i++){
		atom_shared_ptr_vector this_vec = build_peptide_bond(i, pose_);
		rtn_vec.insert(rtn_vec.end(), this_vec.begin(), this_vec.end());
	}

	// even number alphabet build CBs from first_res to last_res
	for (int i = first_res; i <= last_res; i++){
		if (!(pose_->get_residue(i)->get_type() == residue_type("GLY"))){
			if (!(pose_->get_residue(i)->get_type() == residue_type("PRO"))){
				//IC N    C    *CA  CB    1.4592 114.4400  123.2300 111.0900  1.5461
				atom_shared_ptr ca_atm = pose_->get_bb_atom(POSE::CA, i);
				atom_shared_ptr c_atm = pose_->get_bb_atom(POSE::C, i);
				atom_shared_ptr n_atm = pose_->get_bb_atom(POSE::N, i);
				if (ca_atm->isSet() && c_atm->isSet() && n_atm->isSet()){
					const vector3d cb_pos = UTILS::dihedralEnd( n_atm->get_coords(),
							c_atm->get_coords(),
							ca_atm->get_coords(),
							1.5461,
							degrees_to_radians(111.0900),
							degrees_to_radians(123.2300));
					atom_shared_ptr nat = pose_->add_new_atom(cb_pos, atom_type("CB"),i);
					nat->set_b_factor(default_b_factor);
					rtn_vec.push_back(nat);
				}
			}
			else {
				//IC N    C    *CA  CB    1.4585 110.8600  113.7400 111.7400  1.5399
				atom_shared_ptr ca_atm = pose_->get_bb_atom(POSE::CA, i);
				atom_shared_ptr c_atm = pose_->get_bb_atom(POSE::C, i);
				atom_shared_ptr n_atm = pose_->get_bb_atom(POSE::N, i);
				if (ca_atm->isSet() && c_atm->isSet() && n_atm->isSet()){
					const vector3d cb_pos = UTILS::dihedralEnd( n_atm->get_coords(),
							c_atm->get_coords(),
							ca_atm->get_coords(),
							1.5399,
							degrees_to_radians(111.7400),
							degrees_to_radians(113.7400));
					atom_shared_ptr nat = pose_->add_new_atom(cb_pos, atom_type("CB"),i);
					nat->set_b_factor(default_b_factor);
					rtn_vec.push_back(nat);
				}
			}

		}
	}

	for (int i = first_res+1; i <= last_res; i++){
		bool is_def = true;
		const vector3d h_pos = POSE_UTILS::get_ideal_HN(pose_, i, is_def);
		if (is_def){
			atom_shared_ptr nat = pose_->add_new_atom(h_pos, atom_type("H"),i);
			nat->set_b_factor(default_b_factor);
			rtn_vec.push_back(nat);
		}
	}


	return rtn_vec;
}

POSE::atom_shared_ptr_vector alphabet_bb_builder::build_backbone( pose_shared_ptr pose_) const{

	atom_shared_ptr_vector rtn_vec;
	rtn_vec.reserve(10 * pose_->get_residue_count());

	for (int i = 0; i <= pose_->get_residue_count()-letter_length; i++){
		atom_shared_ptr_vector this_vec = this->build_fragment(i, pose_);
		rtn_vec.insert(rtn_vec.end(), this_vec.begin(), this_vec.end());
	}

	if (!is_even){
		for (int i = 0; i < pose_->get_residue_count()-1; i++){
			//IC +N   CA   *C   O     1.3558 116.8400  180.0000 122.5200  1.2297
			atom_shared_ptr ca_atm = pose_->get_bb_atom(POSE::CA, i);
			atom_shared_ptr c_atm = pose_->get_bb_atom(POSE::C, i);
			atom_shared_ptr n_atm = pose_->get_bb_atom(POSE::N, i+1);

			if (ca_atm->isSet() && c_atm->isSet() && n_atm->isSet()){
				const vector3d o_pos = UTILS::dihedralEnd( n_atm->get_coords(),
						ca_atm->get_coords(),
						c_atm->get_coords(),
						1.2297,
						degrees_to_radians(122.5200),
						degrees_to_radians(180.0000));
				atom_shared_ptr nat = pose_->add_new_atom(o_pos, atom_type("O"),i);
				nat->set_b_factor(default_b_factor);
				rtn_vec.push_back(nat);
			}
		}
	}
	else {
		// even number alphabet
		for (int i = 0; i < pose_->get_residue_count(); i++){
			if (!(pose_->get_residue(i)->get_type() == residue_type("GLY"))){
				if (!(pose_->get_residue(i)->get_type() == residue_type("PRO"))){
					//IC N    C    *CA  CB    1.4592 114.4400  123.2300 111.0900  1.5461
					atom_shared_ptr ca_atm = pose_->get_bb_atom(POSE::CA, i);
					atom_shared_ptr c_atm = pose_->get_bb_atom(POSE::C, i);
					atom_shared_ptr n_atm = pose_->get_bb_atom(POSE::N, i);
					if (ca_atm->isSet() && c_atm->isSet() && n_atm->isSet()){
						const vector3d cb_pos = UTILS::dihedralEnd( n_atm->get_coords(),
								c_atm->get_coords(),
								ca_atm->get_coords(),
								1.5461,
								degrees_to_radians(111.0900),
								degrees_to_radians(123.2300));
						atom_shared_ptr nat = pose_->add_new_atom(cb_pos, atom_type("CB"),i);
						nat->set_b_factor(default_b_factor);
						rtn_vec.push_back(nat);
					}
				}
				else {
					//IC N    C    *CA  CB    1.4585 110.8600  113.7400 111.7400  1.5399
					atom_shared_ptr ca_atm = pose_->get_bb_atom(POSE::CA, i);
					atom_shared_ptr c_atm = pose_->get_bb_atom(POSE::C, i);
					atom_shared_ptr n_atm = pose_->get_bb_atom(POSE::N, i);
					if (ca_atm->isSet() && c_atm->isSet() && n_atm->isSet()){
						const vector3d cb_pos = UTILS::dihedralEnd( n_atm->get_coords(),
								c_atm->get_coords(),
								ca_atm->get_coords(),
								1.5399,
								degrees_to_radians(111.7400),
								degrees_to_radians(113.7400));
						atom_shared_ptr nat = pose_->add_new_atom(cb_pos, atom_type("CB"),i);
						nat->set_b_factor(default_b_factor);
						rtn_vec.push_back(nat);
					}
				}

			}
		}
	}

	for (int i = 0; i < pose_->get_residue_count(); i++){
		bool is_def = true;
		const vector3d h_pos = POSE_UTILS::get_ideal_HN(pose_, i, is_def);
		if (is_def){
			atom_shared_ptr nat = pose_->add_new_atom(h_pos, atom_type("H"),i);
			nat->set_b_factor(default_b_factor);
			rtn_vec.push_back(nat);
		}
	}

	//POSE_UTILS::quick_add_HN(pose_, true);
	pose_->index();

	return rtn_vec;

}


}
}
}




