//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * chain.cpp
 *
 *  Created on: Feb 1, 2010
 *      Author: jmacdon
 */

#include "chain.h"
#include "pose.h"


namespace PRODART {
namespace POSE {


boost::shared_ptr<chain> new_chain(){
	boost::shared_ptr<chain> nch(new chain);
	return nch;
}


chain::chain(){

	this->init();

}

chain::chain(const chain& old_chain){

	this->chainID = old_chain.chainID;
	first_internal_residue_index = old_chain.first_internal_residue_index;
	last_internal_residue_index = old_chain.last_internal_residue_index;
	first_bb_atom_index  = old_chain.first_bb_atom_index;
	last_bb_atom_index =  old_chain.last_bb_atom_index;
	b_is_peptide = old_chain.b_is_peptide;


}

chain_shared_ptr chain::clone() const{
	boost::shared_ptr<chain> nch(new chain(*this));
	return nch;
}

void chain::init(){

	first_internal_residue_index = -1;
	last_internal_residue_index = -1;

	first_bb_atom_index = -1;
	last_bb_atom_index = -1;

	chainID = ' ';

	b_is_peptide = true;

}

void chain::set_pose(boost::shared_ptr<pose> ptr){
	this->this_pose = ptr;
}

bool chain::isPeptide() const{
	return this->b_is_peptide;
}

void chain::setIsPeptide(bool val){
	b_is_peptide = val;
}

UTILS::vector3d chain::get_ca_pos(int intindex) const{
	return(this->this_pose.lock()->get_bb_atom(CA, this->first_internal_residue_index+intindex)->get_coords());
}
void chain::set_ca_pos(const UTILS::vector3d pos, int intindex) const{
	this->this_pose.lock()->get_bb_atom(CA, this->first_internal_residue_index+intindex)->set_coords(pos);
}


atom_shared_ptr chain::get_ca(int intindex){
	return(this->this_pose.lock()->get_bb_atom(CA, this->first_internal_residue_index+intindex));
}
const_atom_shared_ptr chain::get_ca(int intindex) const{
	return(this->this_pose.lock()->get_bb_atom(CA, this->first_internal_residue_index+intindex));
}


atom_shared_ptr chain::get_bb_atom(const BBAtomType bb_at_type, int intindex){
	return(this->this_pose.lock()->get_bb_atom(bb_at_type, this->first_internal_residue_index+intindex));
}
const_atom_shared_ptr chain::get_bb_atom(const BBAtomType bb_at_type, int intindex) const{
	return(this->this_pose.lock()->get_bb_atom(bb_at_type, this->first_internal_residue_index+intindex));
}

residue_shared_ptr chain::get_residue(int intindex){
	return(this->this_pose.lock()->get_residue(this->first_internal_residue_index+intindex));
}
const residue_shared_ptr chain::get_residue(int intindex) const{
	return(this->this_pose.lock()->get_residue(this->first_internal_residue_index+intindex));
}


char chain::getChainID() const{
	return this->chainID;
}
void chain::setChainID(const char val){
	this->chainID = val;
}

int chain::length()const{
	return (this->last_internal_residue_index-this->first_internal_residue_index)+1;
}


int chain::get_first_internal_residue_index() const{
	return this->first_internal_residue_index;
}
void chain::set_first_internal_residue_index( const int val){
	this->first_internal_residue_index = val;
}

int chain::get_last_internal_residue_index() const{
	return this->last_internal_residue_index;
}
void chain::set_last_internal_residue_index( const int val){
	this->last_internal_residue_index = val;
}

int chain::get_first_bb_atom_index() const{
	return this->first_bb_atom_index;
}
void chain::set_first_bb_atom_index( const int val){
	this->first_bb_atom_index = val;
}

int chain::get_last_bb_atom_index() const{
	return this->last_bb_atom_index;
}
void chain::set_last_bb_atom_index( const int val){
	this->last_bb_atom_index = val;
}


const residue_type_vector chain::get_sequence() const{
	const_pose_shared_ptr pose_ = this->this_pose.lock();
	residue_type_vector rtn_vec;
	for (int i = this->get_first_internal_residue_index(); i <= this->get_last_internal_residue_index(); i++){
		rtn_vec.push_back(pose_->get_residue(i)->get_type());
	}
	return rtn_vec;
}

const residue_type_vector chain::get_sequence_with_estimated_gaps(const double cutoff) const{
	const_pose_shared_ptr pose_ = this->this_pose.lock();
	residue_type_vector rtn_vec;
	for (int i = this->get_first_internal_residue_index(); i <= this->get_last_internal_residue_index(); i++){
		rtn_vec.push_back(pose_->get_residue(i)->get_type());
		if (i < this->get_last_internal_residue_index()){
			const_atom_shared_ptr this_ca = pose_->get_bb_atom(POSE::CA, i);
			const_atom_shared_ptr next_ca = pose_->get_bb_atom(POSE::CA, i+1);
			if (this_ca && next_ca){
				if (this_ca->isSet() && next_ca->isSet()){
					const double dist = (this_ca->get_coords() - next_ca->get_coords()).mod();
					if (dist > cutoff){
						rtn_vec.push_back(residue_type(""));
					}
				}
			}
		}
	}
	return rtn_vec;
}


}
}
