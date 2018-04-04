//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose.h
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#ifndef POSE_H_
#define POSE_H_

#include <boost/shared_ptr.hpp>
#include "residue.h"
#include "chain.h"
#include "atom.h"
#include "backbone_map.h"
#include "residue_type_map.h"
#include "prodart_env/prodart_env.h"

#include <iostream>
#include <fstream>
#include <boost/math/special_functions/fpclassify.hpp>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>
#include <ctime>
#include <cmath>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
#include<boost/tuple/tuple.hpp>


namespace PRODART {
namespace POSE {

class pose;

typedef boost::shared_ptr<pose> pose_shared_ptr;
typedef boost::shared_ptr<const pose> const_pose_shared_ptr;

typedef boost::weak_ptr<pose> pose_weak_ptr;
typedef boost::weak_ptr<const pose> const_pose_weak_ptr;
//typedef std::map<unsigned long, pose_weak_ptr> ulong_pose_weak_ptr_map;

typedef boost::tuple<atom_shared_ptr, atom_shared_ptr, atom_shared_ptr, atom_shared_ptr> atom_shared_ptr4_tuple;
typedef std::vector<atom_shared_ptr4_tuple> atom_shared_ptr4_tuple_vector;

typedef std::map<POSE::atom_shared_ptr, UTILS::vector3d> atom_shared_ptr_vector3d_map;


class pose {

	//! Factory function to get a pointer to a new pose instance - Do not use any other method to create a new pose object (these should not work anyway)
	friend boost::shared_ptr<pose> new_pose();


public:
	~pose();

	//! make a deep copy of pose object - not thoroughly tested yet!
	pose_shared_ptr clone() const;

	void clear();

	//! output REMARK records separately
	void outputREMARKs(std::ostream& output) const;

	//! output REMARK records separately
	void output_appendlines(std::ostream& output) const;

	//! simple load PDB method - optional check for atomid increment
	std::istream& loadPdb(std::istream& input, const bool check_prev_atomid = true);
	//! simple PDB output method
	void outputPdb(std::ostream& output, bool outputCAonly = false) const;

	int get_chain_count() const;
	//! returns chain by internal chain index (indexed from 0)
	chain_shared_ptr get_chain(const int chain_num);
	const_chain_shared_ptr get_chain(const int chain_num) const;

	//! return chain by chain ID - will return first chain with this ID
	chain_shared_ptr get_chain(const char chainID);
	const_chain_shared_ptr get_chain(const char chainID) const;

	int get_residue_count() const;
	//! get residue using internal residue number (index from 0) - i.e. not the number in PDB
	residue_shared_ptr get_residue(const int residue_num);
	const_residue_shared_ptr get_residue(const int residue_num) const;

	//! get residue using PDB residue number. NOTE: this will be slower than other access methods
	residue_shared_ptr get_residue( std::string pdb_residue_num, const char chainID = ' ', const char iCode = ' ' );
	const_residue_shared_ptr get_residue( std::string pdb_residue_num, const char chainID = ' ', const char iCode = ' '  ) const;

	//! get backbone atom count including non-set and non-active atoms
	const int get_bb_atom_count() const;

	//! get backbone atom - NOTE will return atom ptr even if bb coords have not been set
	atom_shared_ptr get_bb_atom(const BBAtomType,
			const int residue_num);
	//! get backbone atom - NOTE will return atom ptr even if bb coords have not been set
	const_atom_shared_ptr get_bb_atom(const BBAtomType,
			const int residue_num) const;
	//! NOTE will return  vector3d even if coords have not been set
	PRODART::UTILS::vector3d get_bb_atom_coords(const BBAtomType,
			const int residue_num) const;

	atom_shared_ptr get_bb_atom(const int bb_atom_index);
	const_atom_shared_ptr get_bb_atom(const int bb_atom_index) const;
	PRODART::UTILS::vector3d get_bb_atom_coords(const int bb_atom_index) const;

	//! slower method but will also return the first sidechain atoms found with that atom_type - you need to check for NULL pointers (even for backbone atoms as NULL will be returned if the bb atom is not set)
	atom_shared_ptr get_atom(const atom_type type,
			const int residue_num);
	//! slower method but will also return the first sidechain atoms found with that atom_type - you need to check for NULL pointers (even for backbone atoms as NULL will be returned if the bb atom is not set)
	const_atom_shared_ptr get_atom(const atom_type type,
			const int residue_num) const;

	//! returns vector3d coords but will return the value vector3d() if get_atom(residue_num, type) returns NULL so you need to be careful if the atom doesn't exist or are not set
	PRODART::UTILS::vector3d get_atom_coords(const atom_type type,
			const int residue_num) const;

	//! slower method but will also return the first sidechain atoms found with that atom_type - you need to check for NULL pointers (even for backbone atoms as NULL will be returned if the bb atom is not set)
	atom_shared_ptr get_atom(const atom_type type,
			std::string pdb_residue_num, const char chainID = ' ', const char iCode = ' ' );
	//! slower method but will also return the first sidechain atoms found with that atom_type - you need to check for NULL pointers (even for backbone atoms as NULL will be returned if the bb atom is not set)
	const_atom_shared_ptr get_atom(const atom_type type,
			std::string pdb_residue_num, const char chainID = ' ', const char iCode = ' ' ) const;

	PRODART::UTILS::vector3d get_atom_coords(const atom_type type,
			std::string pdb_residue_num, const char chainID = ' ', const char iCode = ' ' ) const;


	void set_bb_atom_coords(const PRODART::UTILS::vector3d& vec,
			const BBAtomType,
			const int residue_num);
	void set_bb_atom_coords(const PRODART::UTILS::vector3d& vec,
			const int bb_atom_index);
	void set_atom_coords(const PRODART::UTILS::vector3d& vec,
			const atom_type type,
			const int residue_num);
	void set_atom_coords(const PRODART::UTILS::vector3d& vec,
			const atom_type type,
			std::string pdb_residue_num, const char chainID = ' ', const char iCode = ' ' );

	atom_shared_ptr add_new_atom(const PRODART::UTILS::vector3d& vec,
			const atom_type type,
			const int residue_num);

	//! Adds n term residue. Seems OK. caution - in testing
	residue_shared_ptr add_nterm_residue(const residue_type type, const int chain_num);
	//! Add c term residue. Seems OK. caution - in testing
	residue_shared_ptr add_cterm_residue(const residue_type type, const int chain_num);

	//! testing don't use. insert residue before internal_residue_num
	residue_shared_ptr insert_residue(const residue_type type, const int internal_residue_num);

	//! make new chain
	chain_shared_ptr add_new_chain(const char chainID, bool peptide_chain = true);
	//! make new chain - auto assign ID
	chain_shared_ptr add_new_chain(bool peptide_chain = true);

	chain_shared_ptr add_duplicated_chain(const_chain_shared_ptr chain_to_copy);

	// TODO merge append_residue() and add_cterm_residue()
	//! append residue to c-term end of chain
	residue_shared_ptr append_residue(const residue_type type, const char chainID);

	//! Deletes residue. Seems OK. caution - testing
	void delete_residue(const unsigned int internal_residue_num);


	//! be careful if pose index is not up to date
	int get_all_atom_count() const;
	//! be careful if pose index is not up to date
	atom_shared_ptr get_atom(const int all_atom_index);
	//! be careful if pose index is not up to date
	const_atom_shared_ptr get_atom(const int all_atom_index) const;

	//! renumber (external) residue numbers in a chain
	void renumber_residues(const int chain_num, const int start_res_num);
	//! renumber (external) residue numbers in a chain
	void renumber_residues(const char chainID, const int start_res_num);
	//! renumber (external) residue numbers for each chain starting from 1
	void renumber_residues(const int start_res_num = 1);

	//! renumber (external) residue numbers in a chain
	void cyclic_renumber_residues(const int chain_num,
			const int start_res_num, const int lower_bound, const int upper_bound, const int period);
	//! renumber (external) residue numbers for each chain starting from 1
	void cyclic_renumber_residues(const int start_res_num, const int lower_bound, const int upper_bound, const int period);

	//! auto renumber chain IDs from A-Z
	void renumber_chainIDs();
	//! auto assign missing chain IDs from A-Z a-z
	void auto_assign_missing_chainIDs();


	atom_shared_ptr4_tuple get_phi_atoms(const int residue_num);
	atom_shared_ptr4_tuple get_psi_atoms(const int residue_num);
	atom_shared_ptr4_tuple get_omega_atoms(const int residue_num);

	//! returns phi dihedral angle in RADIANS - default 0 if required atoms not set
	double get_phi(const int residue_num) const;
	//! returns psi dihedral angle in RADIANS - default 0 if required atoms not set
	double get_psi(const int residue_num) const;
	//! returns omega dihedral angle between this and previous residue in RADIANS - default 0 if required atoms not set
	double get_omega_to_prev(const int residue_num) const;
	//! returns omega dihedral angle between this and next residue in RADIANS - default 0 if required atoms not set
	double get_omega_to_next(const int residue_num) const;



	//! add REMARK record
	void add_remark(const std::string& str);
	//! add REMARK records with "PRODART2" label auto prepended
	void add_remark_with_prodart_label(const std::string& str);
	//! clear REMARK records
	void clear_remarks();
	//! suppress REMARK record output
	void suppress_remark_output();
	//! unsupress REMARK record output
	void unsuppress_remark_output();

	//! add line to append to output
	void add_appendline(std::string str);
	//! clear REMARK records
	void clear_appendlines();
	//! suppress appendline output
	void suppress_appendlines_output();
	//! unsupress appendline output
	void unsuppress_appendlines_output();

	// just write out CAs (toggle) (MIS)
	//bool CAonly;

	//! updates index numbering in the various objects "owned" by _this pose - should be mostly automatically called when necessary
	void index() const;

	//! THIS IS NOT SAFE!!! store a backup of current atom coordinates
	//void backup_coords();
	//! THIS IS NOT SAFE!!! restore backed up atom coordinates
	//void restore_backed_up_coords();

	//! safe replacement for backup_coords()
	atom_shared_ptr_vector3d_map get_all_coords();
	//! safe replacement for restore_backed_up_coords()
	bool restore_coords(atom_shared_ptr_vector3d_map& backup );

	unsigned long get_instance_num() const;
	static unsigned long get_instance_count();

	//static const_pose_shared_ptr get_instance(unsigned long index);

	static int get_num_protein_backbone_atom_types();

	std::string get_label() const;
	void set_label( std::string new_label);


	//! check for NaN and inf in coords
	bool coords_numerically_ok() const;

	void reset_cryst_record();
	void outputCRYST1(std::ostream& output) const;

    int get_model_number() const;
    void set_model_number(int modelNumber);

private:

	void init();
	pose();
	pose(const pose&);



	void outputPdb_ATOM_line(std::ostream& output,
			const atom_shared_ptr currAtom,
			const residue_shared_ptr currResidue ,
			const chain_shared_ptr currChain,
			const int atom_seq_num) const;

	void set_version_date_strings();

	atom_shared_ptr_vector prot_backbone_atom_vec;
	//! mutable index of all atoms - must make sure index is up to date to use functions that rely on this
	mutable atom_shared_ptr_vector all_atoms_vec;
	residue_shared_ptr_vector residue_vec;
	chain_shared_ptr_vector chain_vec;



	const prot_backbone_map*  bb_map;
	const residue_type_map*  rt_map;


	std::string pdbCode;
	std::string label;

	std::vector<std::string> REMARK_vec;
	bool suppress_REMARK_output_flag;

	std::vector<std::string> appendlines_vec;
	bool suppress_appendlines_output_flag;

	std::string version_string;
	std::string run_time_string;

	std::string ASTRAL_version, SCOP_sid, SCOP_sun, SCOP_sccs, Source_PDB, Source_PDB_REVDAT, Region, ASTRAL_SPACI,
		ASTRAL_AEROSPACI, Data_updated_release;

	pose_weak_ptr _this;

	static unsigned long instance_count;

	unsigned long instance_num;

	//static ulong_pose_weak_ptr_map instance_register;


	double cryst_a, cryst_b, cryst_c, cryst_alpha, cryst_beta, cryst_gamma;
	std::string cryst_space_group;
	int cryst_z_value;

	int model_number;


};

//! Factory function to get a pointer to a new pose instance - Do not use any other method to create a new pose object (these should not work anyway)
boost::shared_ptr<pose> new_pose();






/*
 * *****************************************
 * INLINE functions ************************
 * *****************************************
 */


inline residue_shared_ptr pose::get_residue(const int residue_num){
	return residue_vec[residue_num];
}

inline const_residue_shared_ptr pose::get_residue(const int residue_num) const{
	return residue_vec[residue_num];
}

inline atom_shared_ptr pose::get_bb_atom(const BBAtomType bb_at_type,
		const int residue_num){

	return prot_backbone_atom_vec[(residue_num * this->bb_map->get_num_protein_backbone_atoms()) + this->bb_map->get_relative_location(bb_at_type)];

}
inline const_atom_shared_ptr pose::get_bb_atom(const BBAtomType bb_at_type,
		const int residue_num) const{

	return prot_backbone_atom_vec[(residue_num * this->bb_map->get_num_protein_backbone_atoms()) + this->bb_map->get_relative_location(bb_at_type)];

}


inline PRODART::UTILS::vector3d pose::get_bb_atom_coords(const BBAtomType bb_at_type,
		const int residue_num) const{
	return this->get_bb_atom(bb_at_type, residue_num)->get_coords();
}


inline atom_shared_ptr pose::get_bb_atom(const int bb_atom_index){
	return this->prot_backbone_atom_vec[bb_atom_index];
}
inline const_atom_shared_ptr pose::get_bb_atom(const int bb_atom_index) const{
	return this->prot_backbone_atom_vec[bb_atom_index];
}

inline PRODART::UTILS::vector3d pose::get_bb_atom_coords(const int bb_atom_index) const{
	return this->get_bb_atom(bb_atom_index)->get_coords();
}

inline void pose::set_bb_atom_coords(const PRODART::UTILS::vector3d& vec,
		const BBAtomType bb_at_type,
		const int residue_num){
	this->get_bb_atom(bb_at_type, residue_num)->set_coords(vec);
}

inline void pose::set_bb_atom_coords(const PRODART::UTILS::vector3d& vec,
		const int bb_atom_index){
	this->get_bb_atom(bb_atom_index)->set_coords(vec);
}

inline atom_shared_ptr pose::get_atom(const int all_atom_index){
	return this->all_atoms_vec[all_atom_index];
}

inline const_atom_shared_ptr pose::get_atom(const int all_atom_index) const{
	return this->all_atoms_vec[all_atom_index];
}


}
}

#endif /* POSE_H_ */
