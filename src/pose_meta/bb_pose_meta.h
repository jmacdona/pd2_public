/*
 * bb_pose_meta.h
 *
 *  Created on: 22 Oct 2010
 *      Author: jmacdona
 */

#ifndef BB_POSE_META_H_
#define BB_POSE_META_H_

#include <boost/shared_ptr.hpp>
#include "pose_meta_interface.h"
#include "pose_meta_defs.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include "prodart_env/prodart_env.h"
#include <exception>


#include "potentials/restraints_store.h"

namespace PRODART {
namespace POSE {
namespace META {



typedef std::vector<double> double_vector;

class bb_pose_meta;

typedef boost::shared_ptr<bb_pose_meta> bb_pose_meta_shared_ptr;
typedef boost::weak_ptr<bb_pose_meta> bb_pose_meta_weak_ptr;


//! class for storing CA pose meta data (pair lists etc) for use during simulations
class bb_pose_meta : public pose_meta_interface {

	friend bb_pose_meta_shared_ptr new_bb_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach);

private:

	static bool std_atoms_initialised;
	static PRODART::POSE::atom_type bb_aty_N_;
	static PRODART::POSE::atom_type bb_aty_O_;
	static PRODART::POSE::atom_type bb_aty_H_;
	static PRODART::POSE::atom_type bb_aty_CA_;
	static PRODART::POSE::atom_type bb_aty_CB_;
	static PRODART::POSE::atom_type bb_aty_C_;


public:

	//these may not be strictly thread safe but sort out later
	static const PRODART::POSE::atom_type bb_aty_N(){
		return std_atoms_initialised ? bb_aty_N_ : PRODART::POSE::atom_type("N");
	}
	static const PRODART::POSE::atom_type bb_aty_O(){
		return std_atoms_initialised ? bb_aty_O_ : PRODART::POSE::atom_type("O");
	}
	static const PRODART::POSE::atom_type bb_aty_H(){
		return std_atoms_initialised ? bb_aty_H_ : PRODART::POSE::atom_type("H");
	}
	static const PRODART::POSE::atom_type bb_aty_CA(){
		return std_atoms_initialised ? bb_aty_CA_ : PRODART::POSE::atom_type("CA");
	}
	static const PRODART::POSE::atom_type bb_aty_CB(){
		return std_atoms_initialised ? bb_aty_CB_ : PRODART::POSE::atom_type("CB");
	}
	static const PRODART::POSE::atom_type bb_aty_C(){
		return std_atoms_initialised ? bb_aty_C_ : PRODART::POSE::atom_type("C");
	}

protected:

	bb_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach);

	void init();

	void update_pair_lists();

	void initialise_cells();
	void initialise_all_bb_cell();
	void initialise_H_O_cell();

	bool assign_atoms_to_all_bb_sim_cell();
	bool assign_atoms_to_H_O_sim_cell();

	void update_all_bb_pair_list();
	void update_H_O_pair_list();

	void add_sep_info();
	void add_sep_info(PRODART::POSE::META::nb_ele_vector& pair_list);


	void add_bb_bonded_pair(PRODART::POSE::atom_shared_ptr at1,
			PRODART::POSE::atom_shared_ptr at2,
			const double eq_dist,
			const double half_bond_const);
	void createBondedPairList();
	//void updateBondedPairList(bonded_ele_vector& bond_list);

	void add_bb_angle(PRODART::POSE::atom_shared_ptr at1,
			PRODART::POSE::atom_shared_ptr at2,
			PRODART::POSE::atom_shared_ptr at3,
			const double eq_angle,
			const double half_angle_const);
	void createAngleList();
	//void updateAngleList(angle_ele_vector& angle_list);

	void add_bb_dih(dihedral_element& ele);
	void createDihedralList();
	//void updateDihedralList(dihedral_ele_vector& dih_list);


	void createBonded14_15Lists();

	void update_phi_psi_omega();
	void calc_hbond_counts();
	void calc_sec_struct();

	double all_bb_cutoff;
	double all_bb_cell_margin;
	sim_cell all_bb_sim_cell;
	double H_O_cutoff;
	double H_O_cell_margin;
	sim_cell H_O_sim_cell;

	PRODART::POSE::META::nb_ele_vector all_bb_pair_list;
	PRODART::POSE::META::nb_ele_vector H_O_pair_list;

	PRODART::POSE::META::nb_ele_vector bonded_14_list;
	PRODART::POSE::META::nb_ele_vector bonded_15_list;


	double hb_dist_cutoff;

	bonded_ele_vector bb_bond_list;
	angle_ele_vector bb_angle_list;
	dihedral_ele_vector bb_dih_angle_list;

	const prot_backbone_map* const bb_map;

	bb_pose_meta_weak_ptr _this_bb_pose_meta;





public:



	void recalc_pair_lists_dists();

	PRODART::POSE::META::nb_ele_vector& get_all_bb_pair_list();
	PRODART::POSE::META::nb_ele_vector& get_H_O_pair_list();

	PRODART::POSE::META::nb_ele_vector& get_bonded_14_list();
	PRODART::POSE::META::nb_ele_vector& get_bonded_15_list();

	PRODART::POSE::META::bonded_ele_vector& get_bb_bond_list();
	PRODART::POSE::META::angle_ele_vector& get_bb_angle_list();
	PRODART::POSE::META::dihedral_ele_vector& get_bb_dih_angle_list();


	void regularise_atom_names();
	static const double num_bb_atoms;

	poseMetaType get_pose_meta_type() const{
			return pm_bb_pose_meta;
	}


};




bb_pose_meta_shared_ptr new_bb_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach);

}
}
}



#endif /* BB_POSE_META_H_ */
