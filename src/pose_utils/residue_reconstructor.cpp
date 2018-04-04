//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * residue_reconstructor.cpp
 *
 *  Created on: Feb 15, 2010
 *      Author: jmacdona
 */

#include "residue_reconstructor.h"

using namespace PRODART::POSE;
using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::pow;
using std::sqrt;
using std::asin;
using std::string;

using boost::lexical_cast;
using boost::split;
using boost::is_any_of;
using boost::trim;

using std::ios;
using std::setiosflags;
using std::resetiosflags;
using std::setprecision;
using std::setw;

using PRODART::UTILS::vector3d;
using namespace PRODART::POSE;
using namespace PRODART::ROTAMERS;
using namespace std;

typedef std::vector< std::string > string_vector;

namespace PRODART {
namespace POSE_UTILS {


//const prot_backbone_map* const residue_reconstructor::bb_map(prot_backbone_map::Instance());

bool check_bool_vec_all_true(const std::vector<bool> &vec){
	std::vector<bool>::const_iterator iter;

	for (iter = vec.begin(); iter != vec.end(); iter++){
		if (*iter == false){
			return false;
		}
	}

	return true;
}


residue_reconstructor::residue_reconstructor() : bb_map(prot_backbone_map::Instance()){

}

residue_reconstructor::residue_reconstructor(const residue_reconstructor&) : bb_map(prot_backbone_map::Instance()){

}

residue_reconstructor_shared_ptr new_residue_reconstructor(){
	residue_reconstructor_shared_ptr ptr(new residue_reconstructor);
	ptr->_this = ptr;
	return ptr;
}

void residue_reconstructor::make_coor_store(POSE::pose_shared_ptr protein,
			const int internal_residue_num,
			atom_type_temp_store_map& constr_store,
			const bool backbone_only ,
			const bool keep_CB ) const{
	constr_store.clear();
	const vector<BBAtomType> bb_vec = bb_map->get_full_bb_atom_list();

	for (vector<BBAtomType>::const_iterator it = bb_vec.begin(); it != bb_vec.end(); it++){
		atom_shared_ptr atm = protein->get_bb_atom(*it, internal_residue_num);
		if (atm->get_type() == atom_type("CB") && keep_CB){
			constr_store[atm->get_type()] = temp_store(atm->get_coords());
		}
		else if (atm->isSet()){
			constr_store[atm->get_type()] = temp_store(atm->get_coords());
		}
	}
	if (!backbone_only){
		sidechain_shared_ptr sc = protein->get_residue(internal_residue_num)->get_sidechain();
		const int atm_count = sc->get_atom_count();
		for (int i = 0; i  < atm_count; i++){
			atom_shared_ptr atm = sc->get_atom(i);
			if (atm->isSet()){
				constr_store[atm->get_type()] = temp_store(atm->get_coords());
			}
		}
	}
}

void residue_reconstructor::apply_coor_store(POSE::pose_shared_ptr protein,
			const int internal_residue_num,
			atom_type_temp_store_map& constr_store) const{
	const vector<BBAtomType> bb_vec = bb_map->get_full_bb_atom_list();

	for (vector<BBAtomType>::const_iterator it = bb_vec.begin(); it != bb_vec.end(); it++){
		atom_shared_ptr atm = protein->get_bb_atom(*it, internal_residue_num);
		atm->setSet(false);
	}
	sidechain_shared_ptr sc = protein->get_residue(internal_residue_num)->get_sidechain();
	const int atm_count = sc->get_atom_count();
	for (int i = 0; i  < atm_count; i++){
		atom_shared_ptr atm = sc->get_atom(i);
		atm->setSet(false);
	}

	for (atom_type_temp_store_map::const_iterator it = constr_store.begin(); it != constr_store.end(); it++){
		protein->add_new_atom(it->second.coords, it->first, internal_residue_num);
	}


	/*
	const int atm_count = sc->get_atom_count();
	for (int i = 0; i  < atm_count; i++){
		atom_shared_ptr atm = sc->get_atom(i);

	}
	*/
}


void residue_reconstructor::apply_coor_store_to_sc(POSE::pose_shared_ptr protein,
			const int internal_residue_num,
			atom_type_temp_store_map& constr_store,
			ROTAMERS::sidechain_rotamer_shared_ptr sc) const{

	/*
	const vector<BBAtomType> bb_vec = bb_map->get_full_bb_atom_list();

	for (vector<BBAtomType>::const_iterator it = bb_vec.begin(); it != bb_vec.end(); it++){
		atom_shared_ptr atm = protein->get_bb_atom(*it, internal_residue_num);
		atm->setSet(false);
	}
	sidechain_shared_ptr sc = protein->get_residue(internal_residue_num)->get_sidechain();
	const int atm_count = sc->get_atom_count();
	for (int i = 0; i  < atm_count; i++){
		atom_shared_ptr atm = sc->get_atom(i);
		atm->setSet(false);
	}
	*/

	sc->clear();

	std::map<POSE::BBAtomType, boost::tuple<UTILS::vector3d, POSE::atom_type> > changes;
	std::map<POSE::BBAtomType, bool> mask;

	const vector<BBAtomType> bb_vec = bb_map->get_full_bb_atom_list();

	for (vector<BBAtomType>::const_iterator it = bb_vec.begin(); it != bb_vec.end(); it++){
		mask[*it] = false;
	}


	for (atom_type_temp_store_map::const_iterator it = constr_store.begin(); it != constr_store.end(); it++){
		//protein->add_new_atom(it->second.coords, it->first, internal_residue_num);
		if (!bb_map->is_backbone_atom( it->first)){
			sc->add_new_atom(it->second.coords, it->first);
		}
		else {
			const BBAtomType type = bb_map->get_BBAtomType_from_atom_type(it->first);
			mask[type] = true;
			atom_shared_ptr atm = protein->get_bb_atom(type, internal_residue_num);

				//if (atm->isSet() == false || atm->get_type() != it->first || atm->get_coords() != it->second.coords){
			changes[type] = boost::tuple<UTILS::vector3d, POSE::atom_type>(it->second.coords, it->first);
					// is a change
				//}

		}
	}



	sc->set_bacbone_changes(changes);
	sc->set_backbone_mask(mask);
}


int residue_reconstructor::reconstruct_missing_atoms(POSE::pose_shared_ptr protein,
		const int internal_residue_num,
		atom_type_temp_store_map& constr_store) const{
	POSE::residue_shared_ptr res = protein->get_residue(internal_residue_num);

	if (!res) {
		// residue not found!!!
		return false;
	}

	const int num_entries = residue_def.size();
	std::vector<bool> groupOK_vec( num_entries, false );
	int passes = 0;
	bool any_change = false;

	int num_atoms_added = 0;

	while (((check_bool_vec_all_true(groupOK_vec) == false && any_change == true)
			|| passes == 0) && passes < num_entries){

		for (int i = 0; i < num_entries; i++){
			//cout << "entry: " << i << endl;
			bool resAdded = false;
			groupOK_vec[i] = this->add_group(protein,
					internal_residue_num,
					residue_def[i],
					resAdded,
					constr_store);
			if (resAdded == true) num_atoms_added++;
			any_change = any_change || resAdded;
			//cout << "returned:\t" << groupOK_vec[i] << endl;
		}

		passes++;
	}


	//cout << "made num passes " << passes << endl;
	//cout << "changes: " << any_change << endl;
	//cout << "all groupsOK " << check_bool_vec_all_true(groupOK_vec) << endl;
	//cout << "atoms added: " << num_atoms_added << endl;

	protein->index();

	return num_atoms_added;
}

bool residue_reconstructor::reconstruct_missing_atoms(POSE::pose_shared_ptr protein,
		const int internal_residue_num) const{
	POSE::residue_shared_ptr res = protein->get_residue(internal_residue_num);

	if (!res) {
		// residue not found!!!
		return false;
	}

	atom_type_temp_store_map constr_store;
	make_coor_store(protein, internal_residue_num, constr_store);
	const int num_atoms_added = reconstruct_missing_atoms(protein, internal_residue_num, constr_store);
	apply_coor_store(protein, internal_residue_num, constr_store);

	res->set_type(this->res_type);

	if (num_atoms_added > 0){
		return true;
	}

	return false;
}



ROTAMERS::sidechain_rotamer_shared_ptr residue_reconstructor::make_sidechain_rotamer(POSE::pose_shared_ptr protein,
		const int internal_residue_num) const{
	POSE::residue_shared_ptr res = protein->get_residue(internal_residue_num);

	if (!res) {
		// residue not found!!!
		return sidechain_rotamer_shared_ptr();
	}

	sidechain_rotamer_shared_ptr sc = new_sidechain_rotamer();



	atom_type_temp_store_map constr_store;
	make_coor_store(protein, internal_residue_num, constr_store, true, false);
	reconstruct_missing_atoms(protein, internal_residue_num, constr_store);
	apply_coor_store_to_sc(protein, internal_residue_num, constr_store, sc);

	//res->set_type(this->res_type);
	sc->set_residue_type(this->res_type);
	sc->set_residue(res);

	const unsigned int chi_len = chi_defs.size();

	for (unsigned int i = 0 ; i < chi_len; i++){
		sc->set_chi_atoms(i+1, chi_defs[i], chi_fwd[i], chi_bwd[i]);
	}



	return sc;

}

bool residue_reconstructor::old_reconstruct_missing_atoms(POSE::pose_shared_ptr protein,
		const int internal_residue_num) const{

	POSE::residue_shared_ptr res = protein->get_residue(internal_residue_num);

	if (!res) {
		// residue not found!!!
		return false;
	}

	const int num_entries = residue_def.size();
	std::vector<bool> groupOK_vec( num_entries, false );
	int passes = 0;
	bool any_change = false;

	int num_atoms_added = 0;

	while (((check_bool_vec_all_true(groupOK_vec) == false && any_change == true)
			|| passes == 0) && passes < num_entries){

		for (int i = 0; i < num_entries; i++){
			//cout << "entry: " << i << endl;
			bool resAdded = false;
			groupOK_vec[i] = this->add_group(protein,
					internal_residue_num,
					residue_def[i],
					resAdded);
			if (resAdded == true) num_atoms_added++;
			any_change = any_change || resAdded;
			//cout << "returned:\t" << groupOK_vec[i] << endl;
		}

		passes++;
	}


	//cout << "made num passes " << passes << endl;
	//cout << "changes: " << any_change << endl;
	//cout << "all groupsOK " << check_bool_vec_all_true(groupOK_vec) << endl;
	//cout << "atoms added: " << num_atoms_added << endl;

	protein->index();

	return any_change;
}


//! NEED TO ADD CHECK FOR BACKBONE ATOMS THAT ARE NOT "SET"!!!!!
bool residue_reconstructor::add_group(POSE::pose_shared_ptr protein,
		const int internal_residue_num,
		const residue_reconstructor_entry& group,
		bool& res_added) const{

	res_added = false;

	POSE::residue_shared_ptr res = protein->get_residue(internal_residue_num);
	POSE::atom_shared_ptr atoms[4];
	std::vector<bool> isFromOtherRes(4, false);
	std::vector<string> formatted_label(4);
	const residue_shared_ptr prev_res = res->get_prev_residue();
	const residue_shared_ptr next_res = res->get_next_residue();

	/*
	if (!res){
		cout << "res_not_found" << endl;
	}

	if (!prev_res){
		cout << "prev_res_not_found" << endl;
	}

	if (!next_res){
		cout << "next_res_not_found" << endl;
	}
	 */


	for (int i = 0; i < 4; i++){
		string label = group.dihedral_group_atom_names[i];
		trim(label);
		//cout << "\nlabel: " << "'" << label << "'" << endl;
		if (label.substr(0,1).compare("-") == 0){
			// prev residue
			isFromOtherRes[i] = true;
			if (prev_res
					&& internal_residue_num > 0){
				label.erase(0,1);
				atoms[i] = protein->get_atom(atom_type(label),
						internal_residue_num-1);
				//cout << "1.";
			}
			else {
				//cout << "2.";
			}
		}
		else if (label.substr(0,1).compare("+") == 0 ){
			// next residue
			isFromOtherRes[i] = true;
			if (next_res){
				label.erase(0,1);
				atoms[i] = protein->get_atom(atom_type(label),
						internal_residue_num+1);
				//cout << "3.";
			}
			else {
				//cout << "4.";
			}
		}
		else {
			atoms[i] = protein->get_atom(atom_type(label),
					internal_residue_num);
			isFromOtherRes[i] = false;
			//cout << "5.";
		}

		formatted_label[i] = label;
		//cout << endl << "'" << label << "' ";
		if (atoms[i]){
			//cout << "atom_found." << endl;
		}
		else {
			//cout << "atom_not_found." << endl;
		}

	}

	std::vector<bool> atomOK(4, false);

	for (int i = 0; i < 4; i++){
		if (atoms[i]){
			if (atoms[i]->isSet()){
				atomOK[i] = true;
			}
			else {
				atomOK[i] = false;
			}
		}
		else {
			atomOK[i] = false;
		}

	}

	if ((atomOK[0] || isFromOtherRes[0])
	          && (atomOK[1] || isFromOtherRes[1])
	          && (atomOK[2] || isFromOtherRes[2])
	          && (atomOK[3] || isFromOtherRes[3]) ){
		// all atoms already present and correct (set)
		res_added = false;
		return true;
	}
	else if (atomOK[0]
	          && atomOK[1]
	          && atomOK[2]
	          && !atomOK[3]){
		if (isFromOtherRes[3] == true){
			// no need to add missing atoms in other residue
			res_added = false;
			return true;
		}
		// finally we can add missing atom here:

		UTILS::vector3d new_coords;

		//if (!group.isIMPR){
		new_coords = UTILS::dihedralEnd(atoms[0]->get_coords(), atoms[1]->get_coords(), atoms[2]->get_coords(),
				group.bond_3_4,
				group.angle_2_3_4,
				group.dihedral_1_2_3_4);
			/*
			 *
		}
		else {
			new_coords = UTILS::dihedralEnd_d24_a324(atoms[0]->get_coords(), atoms[1]->get_coords(), atoms[2]->get_coords(),
					group.bond_3_4,
					group.angle_2_3_4,
					group.dihedral_1_2_3_4);
		}
		*/

		const int rel_pos = bb_map->get_relative_location(formatted_label[3]);
		//check if is bb atom
		if (rel_pos < bb_map->get_num_protein_backbone_atoms()){
			atom_shared_ptr currAtom = protein->add_new_atom( new_coords,
					atom_type(formatted_label[3]),
					internal_residue_num);
			currAtom->setSet(true);
			currAtom->setActive(true);
			//PRINT_EXPR(atom_type(formatted_label[3]).get_label());
		}
		else {
			atom_shared_ptr currAtom = res->get_sidechain()->add_new_atom(new_coords, atom_type(formatted_label[3]));
			currAtom->setSet(true);
			currAtom->setActive(true);
		}
		res_added = true;
		return true;

	}
	else if (!atomOK[0]
				  && atomOK[1]
	 	          && atomOK[2]
	 	          && atomOK[3]){
		if (isFromOtherRes[0] == true){
			// no need to add missing atoms in other residue
			return true;
		}
		// finally we can add missing atom here:
		UTILS::vector3d new_coords;

		if (!group.isIMPR){
			new_coords = UTILS::dihedralEnd(atoms[3]->get_coords(), atoms[2]->get_coords(), atoms[1]->get_coords(),
					group.bond_1_2,
					group.angle_1_2_3,
					group.dihedral_1_2_3_4);
		}
		else {
			new_coords = UTILS::dihedralEnd_d24_a324(atoms[3]->get_coords(), atoms[2]->get_coords(), atoms[1]->get_coords(),
					group.bond_1_2,					//in fact bond1_3
					group.angle_1_2_3,				// in fact bond_1_3_2
					group.dihedral_1_2_3_4);
		}

		const int rel_pos = bb_map->get_relative_location(formatted_label[0]);
		//check if is bb atom
		if (rel_pos < bb_map->get_num_protein_backbone_atoms()){
			atom_shared_ptr currAtom = protein->add_new_atom( new_coords,
					atom_type(formatted_label[0]),
					internal_residue_num);
			currAtom->setSet(true);
			currAtom->setActive(true);
			//PRINT_EXPR(atom_type(formatted_label[0]).get_label());
		}
		else {
			atom_shared_ptr currAtom = res->get_sidechain()->add_new_atom(new_coords, atom_type());
			currAtom->setSet(true);
			currAtom->setActive(true);
		}
		res_added = true;
		return true;

	}
	else {
		// can't add new atom at this iteration due to lack of atoms
		//cout << "lack of atoms found\n";
		res_added = false;
		return false;
	}

	//cout << "lack of atoms found2\n";
	res_added = false;
	return false;
}


bool residue_reconstructor::add_group(POSE::pose_shared_ptr protein,
		const int internal_residue_num,
		const residue_reconstructor_entry& group,
		bool& res_added,
		atom_type_temp_store_map& constr_store) const{
	res_added = false;

	POSE::residue_shared_ptr res = protein->get_residue(internal_residue_num);
	//POSE::atom_shared_ptr atoms[4];
	UTILS::vector3d coords[4];
	bool atomOK[4] = {false, false, false, false};
	std::vector<bool> isFromOtherRes(4, false);
	std::vector<string> formatted_label(4);
	const residue_shared_ptr prev_res = res->get_prev_residue();
	const residue_shared_ptr next_res = res->get_next_residue();

	/*
	if (!res){
		cout << "res_not_found" << endl;
	}

	if (!prev_res){
		cout << "prev_res_not_found" << endl;
	}

	if (!next_res){
		cout << "next_res_not_found" << endl;
	}
	 */


	for (int i = 0; i < 4; i++){
		string label = group.dihedral_group_atom_names[i];
		trim(label);
		//cout << "\nlabel: " << "'" << label << "'" << endl;
		if (label.substr(0,1).compare("-") == 0){
			// prev residue
			isFromOtherRes[i] = true;
			if (prev_res
					&& internal_residue_num > 0){
				label.erase(0,1);
				atom_shared_ptr atm = protein->get_atom(atom_type(label),
						internal_residue_num-1);
				if (atm){
					if (atm->isSet()){
						atomOK[i] = true;
						coords[i] = atm->get_coords();
					}
				}
				//cout << "1.";
			}
			else {
				atomOK[i] = false;
				//cout << "2.";
			}
		}
		else if (label.substr(0,1).compare("+") == 0 ){
			// next residue
			isFromOtherRes[i] = true;
			if (next_res){
				label.erase(0,1);
				atom_shared_ptr atm = protein->get_atom(atom_type(label),
						internal_residue_num+1);
				if (atm){
					if (atm->isSet()){
						atomOK[i] = true;
						coords[i] = atm->get_coords();
					}
				}
				//cout << "3.";
			}
			else {
				atomOK[i] = false;
				//cout << "4.";
			}
		}
		else {
			/*
			atoms[i] = protein->get_atom(atom_type(label),
					internal_residue_num);
			 */
			isFromOtherRes[i] = false;
			if (constr_store.find(atom_type(label)) != constr_store.end()){
				coords[i] = constr_store[atom_type(label)].coords;
				atomOK[i] = true;
				constr_store[atom_type(label)].is_used_in_def = true;
				//isFromOtherRes[i] = false;
			}
			else {
				atomOK[i] = false;
			}
			//cout << "5.";
		}

		formatted_label[i] = label;


	}

	/*
	std::vector<bool> atomOK(4, false);
	for (int i = 0; i < 4; i++){
		if (atoms[i]){
			if (atoms[i]->isSet()){
				atomOK[i] = true;
			}
			else {
				atomOK[i] = false;
			}
		}
		else {
			atomOK[i] = false;
		}
	}
	*/

	if ((atomOK[0] || isFromOtherRes[0])
	          && (atomOK[1] || isFromOtherRes[1])
	          && (atomOK[2] || isFromOtherRes[2])
	          && (atomOK[3] || isFromOtherRes[3]) ){
		// all atoms already present and correct (set)
		res_added = false;
		return true;
	}
	else if (atomOK[0]
	          && atomOK[1]
	          && atomOK[2]
	          && !atomOK[3]){
		if (isFromOtherRes[3] == true){
			// no need to add missing atoms in other residue
			res_added = false;
			return true;
		}
		// finally we can add missing atom here:

		UTILS::vector3d new_coords;

		//if (!group.isIMPR){
		new_coords = UTILS::dihedralEnd(coords[0], coords[1], coords[2],
				group.bond_3_4,
				group.angle_2_3_4,
				group.dihedral_1_2_3_4);

		/*
		PRINT_EXPR(formatted_label[0]);
		PRINT_EXPR(coords[0]);
		PRINT_EXPR(formatted_label[1]);
		PRINT_EXPR(coords[1]);
		PRINT_EXPR(formatted_label[2]);
		PRINT_EXPR(coords[2]);
		PRINT_EXPR(group.bond_3_4);
		PRINT_EXPR(UTILS::radians_to_degrees(group.angle_2_3_4));
		PRINT_EXPR(UTILS::radians_to_degrees(group.dihedral_1_2_3_4));
		*/

			/*
			 *
		}
		else {
			new_coords = UTILS::dihedralEnd_d24_a324(coords[0], coords[1], coords[2],
					group.bond_3_4,
					group.angle_2_3_4,
					group.dihedral_1_2_3_4);
		}
		*/
		constr_store[atom_type(formatted_label[3])].coords = new_coords;
		constr_store[atom_type(formatted_label[3])].is_used_in_def = true;
		/*
		PRINT_EXPR("**************");
		PRINT_EXPR(formatted_label[3]);
		PRINT_EXPR("**************");
		*/
		/*
		const int rel_pos = bb_map->get_relative_location(formatted_label[3]);
		//check if is bb atom
		if (rel_pos < bb_map->get_num_protein_backbone_atoms()){
			atom_shared_ptr currAtom = protein->add_new_atom( new_coords,
					atom_type(formatted_label[3]),
					internal_residue_num);
			currAtom->setSet(true);
			currAtom->setActive(true);
		}
		else {
			atom_shared_ptr currAtom = res->get_sidechain()->add_new_atom(new_coords, atom_type(formatted_label[3]));
			currAtom->setSet(true);
			currAtom->setActive(true);
		}
		*/
		res_added = true;
		return true;

	}
	else if (!atomOK[0]
				  && atomOK[1]
	 	          && atomOK[2]
	 	          && atomOK[3]){
		if (isFromOtherRes[0] == true){
			// no need to add missing atoms in other residue
			return true;
		}
		// finally we can add missing atom here:
		UTILS::vector3d new_coords;

		if (!group.isIMPR){
			new_coords = UTILS::dihedralEnd(coords[3], coords[2], coords[1],
					group.bond_1_2,
					group.angle_1_2_3,
					group.dihedral_1_2_3_4);
		}
		else {
			new_coords = UTILS::dihedralEnd_d24_a324(coords[3], coords[2], coords[1],
					group.bond_1_2,					//in fact bond1_3
					group.angle_1_2_3,				// in fact bond_1_3_2
					group.dihedral_1_2_3_4);
		}


		constr_store[atom_type(formatted_label[0])].coords = new_coords;
		constr_store[atom_type(formatted_label[0])].is_used_in_def = true;
		//PRINT_EXPR(formatted_label[0]);

		/*
		const int rel_pos = bb_map->get_relative_location(formatted_label[0]);
		//check if is bb atom
		if (rel_pos < bb_map->get_num_protein_backbone_atoms()){
			atom_shared_ptr currAtom = protein->add_new_atom( new_coords,
					atom_type(formatted_label[0]),
					internal_residue_num);
			currAtom->setSet(true);
			currAtom->setActive(true);
		}
		else {
			atom_shared_ptr currAtom = res->get_sidechain()->add_new_atom(new_coords, atom_type());
			currAtom->setSet(true);
			currAtom->setActive(true);
		}
		*/
		res_added = true;
		return true;

	}
	else {
		// can't add new atom at this iteration due to lack of atoms
		//cout << "lack of atoms found\n";
		res_added = false;
		return false;
	}

	//cout << "lack of atoms found2\n";
	res_added = false;
	return false;
}

std::istream& residue_reconstructor::load_residue_def(std::istream& input){

	//file units are degress but we need radians
	const double ang_units = (UTILS::PI / 180.0);

	residue_def.clear();
	chi_defs.clear();
	chi_fwd.clear();
	chi_bwd.clear();


	string lineStr;
	long length, lineNum = 0 ;


	string_vector SplitVec, ProfSplitVec;

	while ( !input.eof() ) {
		getline(input, lineStr);
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length >= 5) {

			split( SplitVec, lineStr, is_any_of("\t ") );
			const string rec_type = SplitVec[0];//lineStr.substr(0,2);

			//cout << endl << lineNum << " " << length << " " << rec_type;

			if (rec_type.compare("IC") == 0 && length >= 66) {

				string an1 = lineStr.substr(3,5);
				string an2 = lineStr.substr(8,5);
				string an3 = lineStr.substr(13,5);
				string an4 = lineStr.substr(18,5);

				trim(an1);
				trim(an2);
				trim(an3);
				trim(an4);


				string b12 = lineStr.substr(23, 7);
				string a123 = lineStr.substr(31, 8);
				string dih1234 = lineStr.substr(40, 9);
				string a234 = lineStr.substr(50, 8);
				string b34 = lineStr.substr(59, 7);

				trim(b12);
				trim(a123);
				trim(dih1234);
				trim(a234);
				trim(b34);

				string isIMPRstr = lineStr.substr(13,1);

				trim(isIMPRstr);

				residue_reconstructor_entry entr;

				if (isIMPRstr.compare("*") == 0){
					entr.isIMPR = true;
					an3.erase(0,1);
				}
				else {
					entr.isIMPR = false;
				}



				entr.dihedral_group_atom_names[0] = an1;
				entr.dihedral_group_atom_names[1] = an2;
				entr.dihedral_group_atom_names[2] = an3;
				entr.dihedral_group_atom_names[3] = an4;

				entr.bond_1_2 = lexical_cast<double>(b12);
				entr.angle_1_2_3 = ang_units * lexical_cast<double>(a123);
				entr.dihedral_1_2_3_4 = ang_units * lexical_cast<double>(dih1234);
				entr.angle_2_3_4 = ang_units * lexical_cast<double>(a234);
				entr.bond_3_4 = lexical_cast<double>(b34);

				residue_def.push_back(entr);

			}
			else if (rec_type.compare("RESI") == 0) {
				//split( SplitVec, lineStr, is_any_of("\t ") );
				const string rstr = SplitVec[1];
				res_type = residue_type(rstr);

			}
			else if (rec_type.compare("CHI_DEF") == 0){
				const unsigned int chi_num = lexical_cast<unsigned int>(SplitVec[1]);
				if (chi_num > chi_defs.size()){
					chi_defs.resize(chi_num);
				}

				atom_type_vector vec;
				for (unsigned int i = 2; i < SplitVec.size(); i++ ){
					string astr = SplitVec[i];
					trim(astr);
					atom_type at(astr);
					vec.push_back(at);
				}
				if (vec.size() != 4){
					cerr << "residue_reconstructor: ERROR: residue definition error in CHI_DEF record" << endl;
				}
				chi_defs[chi_num-1] = vec;
			}
			else if (rec_type.compare("CHI_FWD") == 0){
				const unsigned int chi_num = lexical_cast<unsigned int>(SplitVec[1]);
				if (chi_num > chi_fwd.size()){
					chi_fwd.resize(chi_num);
				}

				atom_type_vector vec;
				for (unsigned int i = 2; i < SplitVec.size(); i++ ){
					string astr = SplitVec[i];
					trim(astr);
					atom_type at(astr);
					vec.push_back(at);
				}
				chi_fwd[chi_num-1] = vec;
			}
			else if (rec_type.compare("CHI_BWD") == 0){
				const unsigned int chi_num = lexical_cast<unsigned int>(SplitVec[1]);
				if (chi_num > chi_bwd.size()){
					chi_bwd.resize(chi_num);
				}

				atom_type_vector vec;
				for (unsigned int i = 2; i < SplitVec.size(); i++ ){
					string astr = SplitVec[i];
					trim(astr);
					atom_type at(astr);
					vec.push_back(at);
				}
				chi_bwd[chi_num-1] = vec;
			}

		}

	}

	if ((chi_defs.size() != chi_fwd.size())
			|| (chi_fwd.size() != chi_bwd.size())){
		cerr << "residue_reconstructor: ERROR: residue definition error in CHI records - each CHI dihedral should have CHI_DEF CHI_FWD CHI_BWD definitions" << endl;
	}

	//cout << residue_def.size() << " definitions loaded\n";

	return input;
}




}
}

