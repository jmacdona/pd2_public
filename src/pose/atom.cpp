//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * atom.cpp
 *
 *  Created on: Feb 1, 2010
 *      Author: jmacdon
 */
#include "atom.h"
#include "residue.h"
#include "sidechain.h"
#include "chain.h"
#include "pose.h"

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

namespace PRODART {
namespace POSE {

boost::shared_ptr<atom> new_atom(){
	boost::shared_ptr<atom> nat(new atom);
	return nat;
}

boost::shared_ptr<atom> new_atom(const atom_type aty, const PRODART::UTILS::vector3d& vec){
	boost::shared_ptr<atom> nat(new atom(aty, vec));
	return nat;
}

//! does not copy pose, residue etc pointers
atom::atom(const atom& old_atom){
	this->coords = old_atom.coords;
	this->type = old_atom.type;
	this->occupancy = old_atom.occupancy;
	this->b_factor = old_atom.b_factor;
	this->active = old_atom.active;
	this->set = old_atom.set;
	this->seq_num = old_atom.seq_num;
	this->charge = old_atom.charge;
	this->mass = old_atom.mass;
	this->has_moved_flag = old_atom.has_moved_flag;
}

//! does not copy pose, residue etc pointers - these need to be arranged afterwards
atom_shared_ptr atom::clone() const{
	boost::shared_ptr<atom> nat(new atom(*this));
	return nat;
}

atom_shared_ptr find( const atom_shared_ptr_vector vec, const atom_type searchtype){

	atom_shared_ptr_vector::const_iterator iter;

	for (iter = vec.begin(); iter != vec.end(); iter++){
		if ((*iter)->get_type() == searchtype){
			return (*iter);
		}
	}

	return atom_shared_ptr();
}


atom::~atom(){
	//cout << "deleting atom: " << type.get_label() << endl;
}

}
}

