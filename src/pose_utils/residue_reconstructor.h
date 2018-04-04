//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * residue_reconstructor.h
 *
 *  Created on: Feb 15, 2010
 *      Author: jmacdona
 */

#ifndef RESIDUE_RECONSTRUCTOR_H_
#define RESIDUE_RECONSTRUCTOR_H_

#include "pose/pose.h"
#include "rotamers/sidechain_rotamer.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>


namespace PRODART {
namespace POSE_UTILS {



class residue_reconstructor;
typedef boost::shared_ptr<residue_reconstructor> residue_reconstructor_shared_ptr;
typedef boost::shared_ptr<const residue_reconstructor> const_residue_reconstructor_shared_ptr;
typedef boost::weak_ptr<residue_reconstructor> residue_reconstructor_weak_ptr;

//! class to recursively reconstruct missing atoms in pose - intending to use charmm parameters
class residue_reconstructor  {

	class residue_reconstructor_entry{

	public:

		std::string dihedral_group_atom_names[4];

		double bond_1_2,  //or bond1_3 when IMPR
			angle_1_2_3,	// or bond_1_3_2 when IMPR
			dihedral_1_2_3_4,
			angle_2_3_4,
			bond_3_4;

		bool isIMPR;



	private:

	};

	typedef std::vector<residue_reconstructor_entry> residue_reconstructor_entry_vector;

	class temp_store{

	public:

		UTILS::vector3d coords;
		bool is_used_in_def;

		temp_store(){
			is_used_in_def = false;
		}
		temp_store(UTILS::vector3d vec){
			coords = vec;
			is_used_in_def = false;
		}

	};

	typedef std::map<POSE::atom_type, temp_store> atom_type_temp_store_map;

	friend residue_reconstructor_shared_ptr new_residue_reconstructor();

protected:

	residue_reconstructor();
	residue_reconstructor(const residue_reconstructor&);

	residue_reconstructor_weak_ptr _this;

	const POSE::prot_backbone_map* const bb_map;
	residue_reconstructor_entry_vector residue_def;
	POSE::residue_type res_type;

	POSE::atom_type_vector_vector chi_defs;
	POSE::atom_type_vector_vector chi_fwd;
	POSE::atom_type_vector_vector chi_bwd;


	void make_coor_store(POSE::pose_shared_ptr protein,
				const int internal_residue_num,
				atom_type_temp_store_map& constr_store,
				const bool backbone_only = false,
				const bool keep_CB = true) const;

	void apply_coor_store(POSE::pose_shared_ptr protein,
				const int internal_residue_num,
				atom_type_temp_store_map& constr_store) const;

	void apply_coor_store_to_sc(POSE::pose_shared_ptr protein,
				const int internal_residue_num,
				atom_type_temp_store_map& constr_store,
				ROTAMERS::sidechain_rotamer_shared_ptr sc) const;

	int reconstruct_missing_atoms(POSE::pose_shared_ptr protein,
			const int internal_residue_num,
			atom_type_temp_store_map& constr_store) const;

	bool old_reconstruct_missing_atoms(POSE::pose_shared_ptr protein,
			const int internal_residue_num) const;

	bool add_group(POSE::pose_shared_ptr protein,
			const int internal_residue_num,
			const residue_reconstructor_entry& group,
			bool& res_added) const;

	bool add_group(POSE::pose_shared_ptr protein,
			const int internal_residue_num,
			const residue_reconstructor_entry& group,
			bool& res_added,
			atom_type_temp_store_map& constr_store) const;

public:


	bool reconstruct_missing_atoms(POSE::pose_shared_ptr protein,
			const int internal_residue_num) const;


	ROTAMERS::sidechain_rotamer_shared_ptr make_sidechain_rotamer(POSE::pose_shared_ptr protein,
			const int internal_residue_num) const;





	//! load internal coordinate definitions (CHARMM-like format) into object
	std::istream& load_residue_def(std::istream& input);

	POSE::residue_type get_residue_type() const{
		return res_type;
	}

};

residue_reconstructor_shared_ptr new_residue_reconstructor();

}
}



#endif /* RESIDUE_RECONSTRUCTOR_H_ */
