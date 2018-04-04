//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose_meta_defs.h
 *
 *  Created on: 4 Jun 2010
 *      Author: jmacdona
 */

#ifndef POSE_META_DEFS_H_
#define POSE_META_DEFS_H_

#include "pose/pose.h"
#include <vector>
#include <list>


namespace PRODART {
namespace POSE {
namespace META {


class coord_rst_element {
public:
	int atom, res_num;
	PRODART::POSE::const_atom_shared_ptr atom_ptr;
	PRODART::POSE::atom_type atype;
	UTILS::vector3d coord;
	double half_bond_const, equilibrium_dist;

	double dist, dist_sq;
	double cached_score;

	coord_rst_element(){
		atom = res_num = 0;
		half_bond_const = equilibrium_dist = cached_score = 0;
		dist = dist_sq = 0;
	}

};

typedef std::vector< coord_rst_element > coord_rst_element_vector;

class pair_element{
public:

	//! this is for contacts
	bool is_valid;

	int atom1 , res_num1;
	int atom2, res_num2;

	PRODART::POSE::const_atom_shared_ptr atom1_ptr;
	PRODART::POSE::const_atom_shared_ptr atom2_ptr;

	PRODART::POSE::atom_type atype1;
	PRODART::POSE::atom_type atype2;

	double cached_score;
	int seq_sep;
	double dist, dist_sq;



	pair_element(){
		dist = dist_sq = 0;
		atom1 = atom2 = 0;
		res_num1 = res_num2 = 0;
		cached_score = 0;
		seq_sep = 999;
		is_valid = true;
	}
};

class nb_pair_element : public pair_element{
	//PRODART::PROT::BondType btype;

public:


	int bond_sep;


	nb_pair_element(){

		bond_sep = 999;
	}

};

class bonded_pair_element : public pair_element {

public:

	double half_bond_const, equilibrium_dist;

	bonded_pair_element() : pair_element(){

		half_bond_const = 0;
		equilibrium_dist = 0;

	}
};

class angle_element {

public:
	int atom1 , res_num1;
	int atom2, res_num2;
	int atom3, res_num3;

	PRODART::POSE::const_atom_shared_ptr atom1_ptr;
	PRODART::POSE::const_atom_shared_ptr atom2_ptr;
	PRODART::POSE::const_atom_shared_ptr atom3_ptr;

	PRODART::POSE::atom_type atype1;
	PRODART::POSE::atom_type atype2;
	PRODART::POSE::atom_type atype3;

	double angle;

	double cached_score;
	double half_angle_const, equilibrium_angle;

	angle_element(){
		atom1 = atom2 = atom3 = 0;
		res_num1 = res_num2 = res_num3 = 0;
		cached_score = 0;
		half_angle_const = 0;
		equilibrium_angle = 0;
		angle = 0;

	}

};

typedef std::vector< angle_element > angle_ele_vector;

class dih_params_element {

public:

	double amplitude;
	double phase;
	double periodicity;
	dih_params_element(){
		amplitude = 0;
		phase = 0;
		periodicity = 0;
	}
	double get_energy(const double dih) const{
		const double energy = 0.5 * this->amplitude * (1+std::cos(this->periodicity*dih - this->phase));
		return energy;
	}
	double get_dE_ddih(const double dih) const{
		const double energy = -0.5 * this->amplitude * this->periodicity * (std::sin(this->periodicity*dih - this->phase));
		return energy;
	}


};
typedef std::vector< dih_params_element > dih_params_ele_vector;


class simple_harmonic_dihedral_element {
public:

	int atom1 , res_num1;
	int atom2, res_num2;
	int atom3, res_num3;
	int atom4, res_num4;

	PRODART::POSE::const_atom_shared_ptr atom1_ptr;
	PRODART::POSE::const_atom_shared_ptr atom2_ptr;
	PRODART::POSE::const_atom_shared_ptr atom3_ptr;
	PRODART::POSE::const_atom_shared_ptr atom4_ptr;

	double dih_angle;

	double equil_dih_angle;
	double weight;
	simple_harmonic_dihedral_element(){
		dih_angle = equil_dih_angle = weight = 0;
		atom1 = atom2 = atom3 = 0;
		res_num1 = res_num2 = res_num3 = 0;
	}

};

typedef std::vector< simple_harmonic_dihedral_element > simple_harmonic_dihedral_element_vector;

class dihedral_element {

public:
	int atom1 , res_num1;
	int atom2, res_num2;
	int atom3, res_num3;
	int atom4, res_num4;

	PRODART::POSE::const_atom_shared_ptr atom1_ptr;
	PRODART::POSE::const_atom_shared_ptr atom2_ptr;
	PRODART::POSE::const_atom_shared_ptr atom3_ptr;
	PRODART::POSE::const_atom_shared_ptr atom4_ptr;

	PRODART::POSE::atom_type atype1;
	PRODART::POSE::atom_type atype2;
	PRODART::POSE::atom_type atype3;
	PRODART::POSE::atom_type atype4;

	double dih_angle;

	double cached_score;
	dih_params_ele_vector params;


	dihedral_element(){
		params.clear();
		params.reserve(10);
		atom1 = atom2 = atom3 = 0;
		res_num1 = res_num2 = res_num3 = 0;
		cached_score = 0;
		dih_angle = 0;

	}

};

typedef std::vector< dihedral_element > dihedral_ele_vector;

template <class T>
class template_pair {
public:
	T first;
	T second;

	template_pair(){
		this->first = 0;
		this->second = 0;
	}

	template_pair(const T val1, const T val2){
		this->first = val1;
		this->second = val2;
	}

	~template_pair(){

	}

	bool operator<( const template_pair<T> & other ) const{
		if (this->first < other.first ) {
			return true;
		}
		else if (this->first == other.first
				&& this->second < other.second){
			return true;
		}

		return false;
	}

	bool operator==( const template_pair<T> & other ) const{
		if (this->first == other.first
				&& this->second == other.second){
			return true;
		}

		return false;
	}

};


typedef std::list< nb_pair_element > nb_ele_list;
typedef std::vector< nb_pair_element > nb_ele_vector;

typedef std::vector< bonded_pair_element > bonded_ele_vector;

typedef std::vector<int> int_vector;

typedef template_pair< int > int_pair;
typedef template_pair< double > double_pair;
typedef std::vector<int_pair> int_pair_vector;

enum fragmentPositionType{fragMIDDLE = 0,
	fragNTERM = 1,
	fragCTERM = 2,
	fragNandCTERM = 3};

class frag4_element{

public:
	PRODART::POSE::const_atom_shared_ptr CA_atoms[4];
	fragmentPositionType frag_pos;
	int frag_type_num;
	PRODART::POSE::four_state_sec_struct frag_ss_class;

	double tau, omega1, omega2;

	int CA_1_residue_number;

	int IBBB_key;

};

typedef std::vector<frag4_element> frag4_vector;


class ca_hbond_collection{

public:

	PRODART::POSE::const_atom_shared_ptr N_i, O_i, CA_i, CA_i_p1;

};
typedef boost::shared_ptr<ca_hbond_collection> ca_hbond_collection_shared_ptr;
typedef std::map<PRODART::POSE::const_atom_shared_ptr, ca_hbond_collection_shared_ptr > const_atom_shared_ptr_ca_hbond_collection_shared_ptr_map;



class sse_axis_element {
public:
	coord_rst_element start, end;
	double seg_intersect_dist_half_bond_const;
	double seg_intersect_dih_half_bond_const;
	double seg_intersect_ang_half_bond_const;
	POSE::four_state_sec_struct secs;

	sse_axis_element(){
		seg_intersect_dist_half_bond_const = 0;
		seg_intersect_dih_half_bond_const = 0;
		seg_intersect_ang_half_bond_const = 0;
		start.half_bond_const = 0;
		end.half_bond_const = 0;
		secs = POSE::ss4_UNDEF;
	}
};

typedef std::vector<sse_axis_element> sse_axis_element_vector;

class upper_lower_bonded_pair_element : public pair_element {

public:

	double half_bond_const,
		equilibrium_dist_lower,
		equilibrium_dist_upper;

	upper_lower_bonded_pair_element() : pair_element(){

		half_bond_const = 0;
		equilibrium_dist_lower = 0;
		equilibrium_dist_upper = 0;

	}
};
typedef std::vector<upper_lower_bonded_pair_element> upper_lower_bonded_pair_element_vector;

}
}
}


#endif /* POSE_META_DEFS_H_ */
