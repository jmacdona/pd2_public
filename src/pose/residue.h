//
// (c)  JAMES T. MACDONALD 2010 and 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * residue.h
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#ifndef RESIDUE_H_
#define RESIDUE_H_
#include "atom.h"
#include "residue_type.h"
#include "sidechain.h"
#include "backbone_map.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include "string"
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
//#include <boost/make_shared.hpp>
#include <cassert>
#include <exception>





namespace PRODART {
namespace POSE {

class pose;
class chain;
class residue;

//! 4 state secondary structure enumeration
enum four_state_sec_struct{
	ss4_OTHER = 0,
	ss4_HELIX = 1,
	ss4_STRAND = 2,
	ss4_CIS = 3,
	ss4_UNDEF
};

//! 4 state secondary structure vector
typedef std::vector<four_state_sec_struct> four_state_sec_struct_vector;

//! 3 state secondary structure enumeration
enum three_state_sec_struct{
	ss3_OTHER = 0,
	ss3_HELIX = 1,
	ss3_STRAND = 2,
	ss3_UNDEF
};

//! 3 state secondary structure vector
typedef std::vector<three_state_sec_struct> three_state_sec_struct_vector;

char to_char(four_state_sec_struct);
char to_char(three_state_sec_struct);

four_state_sec_struct to_4state(const char c);
three_state_sec_struct to_3state(const char c);

four_state_sec_struct to_4state(three_state_sec_struct val);
three_state_sec_struct to_3state(four_state_sec_struct val);

std::ostream &operator<<( std::ostream &output, const four_state_sec_struct_vector& vec);
std::ostream &operator<<( std::ostream &output, const three_state_sec_struct_vector& vec);

std::istream& operator>>( std::istream& input, four_state_sec_struct_vector& vec );
std::istream& operator>>( std::istream& input, three_state_sec_struct_vector& vec );

//enum AtomType {UNDEF_ATOM, N, H, CA, 1HA, 2HA, CB, 1HB, 2HB, 3HB, C, O, END_ATOMTYPE_ENUM};

typedef boost::shared_ptr<residue> residue_shared_ptr;
typedef boost::shared_ptr<const residue> const_residue_shared_ptr;
typedef std::vector<residue_shared_ptr> residue_shared_ptr_vector;

typedef boost::weak_ptr<residue> residue_weak_ptr;
typedef boost::weak_ptr<const residue> const_residue_weak_ptr;


class residue {

	friend boost::shared_ptr<residue> new_residue();
	friend class pose;
	friend class chain;


public:

	boost::shared_ptr<pose> get_pose();
	boost::shared_ptr<const pose> get_pose() const;

	boost::shared_ptr<chain> get_chain();
	boost::shared_ptr<const chain> get_chain() const;

	sidechain_shared_ptr get_sidechain();
	const_sidechain_shared_ptr get_sidechain() const;

	//! replace sidechain
	void swap_sidechain(sidechain_shared_ptr new_sc);

	//! must be used carefully as contained atoms' pointers to pose, residue etc must already be set
	void set_sidechain(sidechain_shared_ptr ptr);

	residue_type get_type() const;
	void set_type(residue_type);

	//! get Code for insertion of residues - sometimes used to show where residues have been inserted in PDB files
	char get_icode() const;
	void set_icode(const char);

	int get_internal_residue_index() const;

	std::string get_pdb_residue_index() const;
	std::string get_trimmed_pdb_residue_index() const;
	void set_pdb_residue_index(std::string str);
	void set_pdb_residue_index(const int num);

	int get_first_bb_atom_index() const;

	residue_shared_ptr get_prev_residue();
	const_residue_shared_ptr get_prev_residue() const;

	residue_shared_ptr get_next_residue();
	const_residue_shared_ptr get_next_residue() const;

	//! is residue an N- or C-terminal residue (based in whether there are links to to previous and next residues)
	bool is_terminal() const;

	//! get backbone atom - NOTE will return atom ptr even if bb coords have not been set
	atom_shared_ptr get_bb_atom(const BBAtomType);
	//! get backbone atom - NOTE will return atom ptr even if bb coords have not been set
	const_atom_shared_ptr get_bb_atom(const BBAtomType) const;

	//! get backbone atom - NOTE will return atom ptr even if bb coords have not been set
	atom_shared_ptr get_atom(const atom_type);
	//! get backbone atom - NOTE will return atom ptr even if bb coords have not been set
	const_atom_shared_ptr get_atom(const atom_type) const;

	//! get vector of all atoms in residue (backbone and sidechain)
	atom_shared_ptr_vector get_all_atoms();
	const_atom_shared_ptr_vector get_all_atoms() const;

	//! returns phi dihedral angle in RADIANS - default 0 if required atoms not set
	double get_phi() const;
	//! returns psi dihedral angle in RADIANS - default 0 if required atoms not set
	double get_psi() const;
	//! returns omega dihedral angle between this and previous residue in RADIANS - default 0 if required atoms not set
	double get_omega_to_prev() const;
	//! returns omega dihedral angle between this and next residue in RADIANS - default 0 if required atoms not set
	double get_omega_to_next() const;

protected:

	residue_shared_ptr clone() const;

	void set_internal_residue_index(const int);
	void set_first_bb_atom_index(const int);

	void set_pose(boost::shared_ptr<pose> ptr);
	void set_chain(boost::shared_ptr<chain> ptr);

	void set_prev_residue(residue_shared_ptr ptr);
	void set_next_residue(residue_shared_ptr ptr);


private:

	residue();
	residue(const residue&);


	void init();

	const prot_backbone_map* const bb_map;

	boost::weak_ptr<pose> this_pose;
	boost::weak_ptr<chain> this_chain;

	residue_weak_ptr _this;
	// previous (N) residue in chain
	residue_weak_ptr prev_residue;
	// next (C) residue in chain
	residue_weak_ptr next_residue;

	int first_bb_atom_index;

	int internal_residue_index;
	//! residue index as read from PDB file
	char pdb_residue_index[4];
	//! insertion code as read from PDB file
	char icode;
	sidechain_shared_ptr sc;

	residue_type type;




};

boost::shared_ptr<residue> new_residue();


inline residue_type residue::get_type() const{
	return type;
}
inline void residue::set_type(residue_type new_type){
	type = new_type;
}

inline int residue::get_first_bb_atom_index() const{
	return first_bb_atom_index;
}

inline boost::shared_ptr<pose> residue::get_pose(){
	return this->this_pose.lock();
}

inline boost::shared_ptr<const pose> residue::get_pose() const{
	return this->this_pose.lock();
}

inline boost::shared_ptr<chain> residue::get_chain(){
	return this->this_chain.lock();
}
inline boost::shared_ptr<const chain> residue::get_chain() const{
	return this->this_chain.lock();
}

inline sidechain_shared_ptr residue::get_sidechain(){
	return this->sc;
}

inline const_sidechain_shared_ptr residue::get_sidechain() const{
	return this->sc;
}

inline residue_shared_ptr residue::get_prev_residue(){
	return this->prev_residue.lock();
}
inline const_residue_shared_ptr residue::get_prev_residue() const{
	return this->prev_residue.lock();
}

inline residue_shared_ptr residue::get_next_residue(){
	return this->next_residue.lock();
}
inline const_residue_shared_ptr residue::get_next_residue() const{
	return this->next_residue.lock();
}

inline bool residue::is_terminal() const{
	if (this->get_next_residue() && this->get_prev_residue()){
		return false;
	}
	else {
		return true;
	}
}

inline void residue::set_sidechain(sidechain_shared_ptr ptr){
	this->sc = ptr;
	sc->set_residue(_this.lock());
}



}
}

#endif /* RESIDUE_H_ */
