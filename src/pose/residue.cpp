//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * residue.cpp
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#include "residue.h"
#include "pose.h"
#include "chain.h"


using boost::trim;
using std::string;

namespace PRODART {
namespace POSE {


//const prot_backbone_map* const residue::bb_map(prot_backbone_map::Instance());

boost::shared_ptr<residue> new_residue(){
	boost::shared_ptr<residue> nres(new residue);
	nres->_this = nres;
	nres->init();
	return nres;
}

residue::residue(const residue& old_residue) : bb_map(prot_backbone_map::Instance()){
	//this->init();
	this->first_bb_atom_index = old_residue.first_bb_atom_index;
	this->internal_residue_index = old_residue.internal_residue_index;
	for (int i = 0; i < 4; i++){
		pdb_residue_index[i] = old_residue.pdb_residue_index[i];
	}
	icode = old_residue.icode;
	type = old_residue.type;

	//copy sidechain
	this->sc = old_residue.sc->clone();

}

//! TODO NOTE: need to sort out sidechain cloning
residue_shared_ptr residue::clone() const{
	boost::shared_ptr<residue> nres(new residue(*this));
	nres->_this = nres;
	nres->init();
	return nres;
}

void residue::init(){
	this->sc = new_sidechain();
	this->sc->set_residue(_this.lock());

}

residue::residue() : bb_map(prot_backbone_map::Instance()){
	first_bb_atom_index = 0;
	internal_residue_index = 0;
	for (int i = 0; i < 4; i++){
		pdb_residue_index[i] = ' ';
	}
	icode = ' ';
	//owner = NULL;
	//this->init();
}

void residue::set_pose(boost::shared_ptr<pose> ptr){
	this->this_pose = ptr;
}

void residue::set_chain(boost::shared_ptr<chain> ptr){
	this->this_chain = ptr;
}


char residue::get_icode() const{
	return this->icode;
}
void residue::set_icode(const char val){
	this->icode = val;
}


int residue::get_internal_residue_index() const{
	return this->internal_residue_index;
}
void residue::set_internal_residue_index(const int val){
	this->internal_residue_index = val;
}

void residue::set_first_bb_atom_index(const int val){
	this->first_bb_atom_index = val;
}

std::string residue::get_pdb_residue_index() const{
	std::string ret_str(4, ' ');
	for (int i = 0; i < 4; i++){
		ret_str[i] = this->pdb_residue_index[i];
	}
	return ret_str;
}

std::string residue::get_trimmed_pdb_residue_index() const{
	std::string ret_str = get_pdb_residue_index();
	trim(ret_str);
	return ret_str;
}

void residue::set_pdb_residue_index(std::string str){
	trim(str);
	const int str_len = str.length();

	if (str_len < 4){
		str.insert(0, 4 - str_len, ' ');
	}
	else if (str_len > 4){
		str.erase(0, str_len - 4);
	}

	for (int i = 3; i >= 0; i--){
		//if (i < str_len){
		this->pdb_residue_index[i] = str[i];
			/*
		}
		else {
			this->pdb_residue_index[i] = ' ';
		}
		*/
	}
}

void residue::set_pdb_residue_index(const int num){
	string str(boost::lexical_cast<string>(num));
	this->set_pdb_residue_index(str);
}

void residue::set_prev_residue(residue_shared_ptr ptr){
	this->prev_residue = ptr;
}

void residue::set_next_residue(residue_shared_ptr ptr){
	this->next_residue = ptr;
}

atom_shared_ptr residue::get_bb_atom(const BBAtomType at){
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_bb_atom(at, res_num);
}

const_atom_shared_ptr residue::get_bb_atom(const BBAtomType at) const{
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_bb_atom(at, res_num);
}

atom_shared_ptr residue::get_atom(const atom_type at){
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_atom(at, res_num);
}
const_atom_shared_ptr residue::get_atom(const atom_type at) const{
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_atom(at, res_num);
}

//! replace sidechain
void residue::swap_sidechain(sidechain_shared_ptr new_sc){
	assert(new_sc);
	this->sc = new_sc;
	assert(sc->get_residue() == _this.lock());
}


atom_shared_ptr_vector residue::get_all_atoms(){
	atom_shared_ptr_vector rtn_vec;
	pose_shared_ptr pose_ =  this->get_pose();
	rtn_vec.reserve(bb_map->get_num_protein_backbone_atoms() + sc->get_atom_count());
	for (int i = this->get_first_bb_atom_index(); i < (this->get_first_bb_atom_index() + bb_map->get_num_protein_backbone_atoms()); i++){
		atom_shared_ptr ptr = pose_->get_bb_atom(i);
		if (ptr->get_residue() != _this.lock()){
			std::cerr << "ERROR: bug in function atom::get_all_atoms() - not getting correct atoms";
		}
		rtn_vec.push_back(ptr);
	}

	atom_shared_ptr_vector sc_rtn_vec = sc->get_all_sc_atoms();

	rtn_vec.insert(rtn_vec.end(), sc_rtn_vec.begin(), sc_rtn_vec.end());

	return rtn_vec;
}

const_atom_shared_ptr_vector residue::get_all_atoms() const{
	const_atom_shared_ptr_vector rtn_vec;
	const_pose_shared_ptr pose_ =  this->get_pose();
	rtn_vec.reserve(bb_map->get_num_protein_backbone_atoms() + sc->get_atom_count());
	for (int i = this->get_first_bb_atom_index(); i < (this->get_first_bb_atom_index() + bb_map->get_num_protein_backbone_atoms()); i++){
		const_atom_shared_ptr ptr = pose_->get_bb_atom(i);
		if (ptr->get_residue() != _this.lock()){
			std::cerr << "ERROR: bug in function atom::get_all_atoms() const - not getting correct atoms";
		}
		rtn_vec.push_back(ptr);
	}

	const_atom_shared_ptr_vector sc_rtn_vec = get_sidechain()->get_all_sc_atoms();

	rtn_vec.insert(rtn_vec.end(), sc_rtn_vec.begin(), sc_rtn_vec.end());

	return rtn_vec;
}


double residue::get_phi() const{
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_phi(res_num);
}

double residue::get_psi() const{
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_psi(res_num);
}

double residue::get_omega_to_prev() const{
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_omega_to_prev(res_num);
}

double residue::get_omega_to_next() const{
	const int res_num = this->get_internal_residue_index();
	return this->get_pose()->get_omega_to_next(res_num);
}

char to_char(four_state_sec_struct val){
	if (val == PRODART::POSE::ss4_OTHER){
		return '-';
	}
	else if (val == PRODART::POSE::ss4_HELIX){
		return 'H';
	}
	else if (val == PRODART::POSE::ss4_STRAND){
		return 'E';
	}
	else if (val == PRODART::POSE::ss4_CIS){
		return 'C';
	}
	return '?';
}

char to_char(three_state_sec_struct val){
	if (val == PRODART::POSE::ss3_OTHER){
		return '-';
	}
	else if (val == PRODART::POSE::ss3_HELIX){
		return 'H';
	}
	else if (val == PRODART::POSE::ss3_STRAND){
		return 'E';
	}
	return '?';
}

four_state_sec_struct to_4state(const char c){
	if (c == '-'){
		return PRODART::POSE::ss4_OTHER;//'-';
	}
	else if (c == 'H'){
		return PRODART::POSE::ss4_HELIX;//'H';
	}
	else if (c == 'E'){
		return PRODART::POSE::ss4_STRAND;//'E';
	}
	else if (c == 'C' ){
		return PRODART::POSE::ss4_CIS;//'C';
	}
	return PRODART::POSE::ss4_UNDEF;
}
three_state_sec_struct to_3state(const char c){
	if (c == '-'){
		return PRODART::POSE::ss3_OTHER;//'-';
	}
	else if (c == 'H'){
		return PRODART::POSE::ss3_HELIX;//'H';
	}
	else if (c == 'E'){
		return PRODART::POSE::ss3_STRAND;//'E';
	}

	return PRODART::POSE::ss3_UNDEF;
}

four_state_sec_struct to_4state(three_state_sec_struct val){
	if (val == PRODART::POSE::ss3_OTHER){
		return PRODART::POSE::ss4_OTHER;//'-';
	}
	else if (val == PRODART::POSE::ss3_HELIX){
		return PRODART::POSE::ss4_HELIX;//'H';
	}
	else if (val == PRODART::POSE::ss3_STRAND){
		return PRODART::POSE::ss4_STRAND;//'E';
	}

	return PRODART::POSE::ss4_UNDEF;
}

three_state_sec_struct to_3state(four_state_sec_struct val){
	if (val == PRODART::POSE::ss4_OTHER){
		return PRODART::POSE::ss3_OTHER;//'-';
	}
	else if (val == PRODART::POSE::ss4_HELIX){
		return PRODART::POSE::ss3_HELIX;//'H';
	}
	else if (val == PRODART::POSE::ss4_STRAND){
		return PRODART::POSE::ss3_STRAND;//'E';
	}
	else if (val == PRODART::POSE::ss4_CIS ){
		return PRODART::POSE::ss3_OTHER;//'C';
	}
	return PRODART::POSE::ss3_UNDEF;
}

class read_secs_exception: public std::exception{
  virtual const char* what() const throw()
  {
    return "read_secs_exception";
  }
};

std::istream& operator>>( std::istream& input, four_state_sec_struct_vector& vec ){
	string lineStr;
	bool isOK = true;
	while ( !input.eof() ) {
		getline(input, lineStr);
		trim(lineStr);

		const long length = lineStr.length();
		for (int i =0; i < length; i++){
			four_state_sec_struct sec = to_4state(lineStr[i]);
			vec.push_back(sec);
			if (sec == ss4_UNDEF){
				isOK = false;
			}
		}
	}

	if (!isOK){
		std::cerr << "ERROR reading four_state_sec_struct_vector: unrecognised characters" << std::endl;
		throw read_secs_exception();
	}

	return input;
}
std::istream& operator>>( std::istream& input, three_state_sec_struct_vector& vec ){
	string lineStr;
	bool isOK = true;
	while ( !input.eof() ) {
		getline(input, lineStr);
		trim(lineStr);

		const long length = lineStr.length();
		for (int i =0; i < length; i++){
			three_state_sec_struct sec = to_3state(lineStr[i]);
			vec.push_back(sec);
			if (sec == ss3_UNDEF){
				isOK = false;
			}
		}
	}

	if (!isOK){
		std::cerr << "ERROR reading three_state_sec_struct_vector: unrecognised characters" << std::endl;
		throw read_secs_exception();
	}

	return input;
}

std::ostream &operator<<( std::ostream &output, const four_state_sec_struct_vector& vec ){
	four_state_sec_struct_vector::const_iterator iter;

	for (iter = vec.begin(); iter != vec.end(); iter++){
		output << to_char(*iter);
	}
	return output;
}

std::ostream &operator<<( std::ostream &output, const three_state_sec_struct_vector& vec ){
	three_state_sec_struct_vector::const_iterator iter;

	for (iter = vec.begin(); iter != vec.end(); iter++){
		output << to_char(*iter);
	}
	return output;
}


}
}
