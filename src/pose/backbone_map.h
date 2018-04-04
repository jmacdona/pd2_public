//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * backbone_map.h
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#ifndef BACKBONE_MAP_H_
#define BACKBONE_MAP_H_

#include "atom_type.h"
#include <boost/shared_ptr.hpp>
#include <map>
#include <vector>
#include <iostream>


#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <boost/thread.hpp>
#include "prodart_env/prodart_env.h"

namespace PRODART {
namespace POSE {


typedef std::map<atom_type, int> atom_type_int_map;
//typedef std::map<int, BBAtomType> int_BBAtomType_map;
typedef std::vector<double> double_vector;
typedef std::vector<int> int_vector;

const int num_protein_backbone_atoms = 14;//11;

/*enum AtomType {
 * N,	0
 * H, 	1
 * CA, 	2
 * 1HA, 3
 * 2HA, 4
 * CB, 	5
 * 1HB, 6
 * 2HB, 7
 * 3HB, 8
 * C, 	9
 * O, 	10
 * };
 */

//! backbone atom types for quick access
enum BBAtomType { N, H, CA, HA1, HA2, CB, HB1, HB2, HB3, C,  O, HA, H1, H2, H3, OXT, END_BBAtomType};



//! singleton class
class prot_backbone_map{

private:
	prot_backbone_map(const prot_backbone_map&);

	static void Init();

	void init_bond_sep_data();

	atom_type_int_map at_map;

	std::vector<int> at_enum_vec;
	std::vector<int> inv_at_enum_vec;


	std::vector<int> pdb_output_order_vec;

	std::vector<BBAtomType> smaller_bb_atom_list;
	std::vector<BBAtomType> small_bb_atom_list;
	std::vector<BBAtomType> full_bb_atom_list;
	std::vector<BBAtomType> bb_hydrogen_atom_list;


	int** internal_residue_bond_sep_mat;
	int_vector bonds_to_N, bonds_to_C;


protected:
	prot_backbone_map();
	virtual ~prot_backbone_map();

public:

	static const prot_backbone_map* Instance();

	static int get_num_protein_backbone_atoms() {
		return num_protein_backbone_atoms;
	}

	bool is_backbone_atom(const atom_type) const;

	int get_relative_location(const atom_type) const;
	int get_relative_location(const BBAtomType) const;
	BBAtomType get_BBAtomType_from_rel_loc(const int relative_location) const;
	BBAtomType get_BBAtomType_from_atom_type(const atom_type) const;

	int get_output_order(const int index) const;

	int get_bond_sep(const int seq_sep, const BBAtomType lower_at, const BBAtomType upper_at) const;
	int get_bond_sep(const int seq_sep, const atom_type lower_at, const atom_type upper_at) const;


	//! returns {N, CA, C, O, CB}
	const std::vector<BBAtomType>& get_smaller_bb_atom_list() const {
		return smaller_bb_atom_list;
	}

	//! returns {N, CA, C, O, CB, H}
	const std::vector<BBAtomType>& get_small_bb_atom_list() const {
		return small_bb_atom_list;
	}

	//! returns full list of bb atoms types
	const std::vector<BBAtomType>& get_full_bb_atom_list() const {
		return full_bb_atom_list;
	}

	//! returns full list of bb atoms types
	const std::vector<BBAtomType>& get_bb_hydrogen_atom_list() const {
		return bb_hydrogen_atom_list;
	}


};

inline BBAtomType prot_backbone_map::get_BBAtomType_from_rel_loc(const int relative_location) const{
	if (relative_location >= num_protein_backbone_atoms || relative_location < 0){
		return END_BBAtomType;
	}

	return static_cast<BBAtomType>(inv_at_enum_vec[relative_location]);

}

inline BBAtomType prot_backbone_map::get_BBAtomType_from_atom_type(const atom_type type) const{
	const int relative_location = this->get_relative_location(type);
	return get_BBAtomType_from_rel_loc(relative_location);
}


inline bool prot_backbone_map::is_backbone_atom(const atom_type type) const{
	atom_type_int_map::const_iterator it = at_map.find(type);
	if (it != at_map.end()){
		return true;
	}
	else {
		return false;
	}
}

inline int prot_backbone_map::get_relative_location(const atom_type type) const{
	atom_type_int_map::const_iterator it = at_map.find(type);
	if (it != at_map.end()){
		return it->second;
	}
	else {
		return num_protein_backbone_atoms;
	}
}

inline int prot_backbone_map::get_relative_location(const BBAtomType at) const{
	return at_enum_vec[static_cast<int>(at)];
}

}
}


#endif /* BACKBONE_MAP_H_ */
