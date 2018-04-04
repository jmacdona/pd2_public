//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * ca_pose_meta.h
 *
 *  Created on: 5 Mar 2010
 *      Author: jmacdona
 */

#ifndef CA_POSE_META_H_
#define CA_POSE_META_H_
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

namespace FRAG_CLASS {
class fragment_classifier_interface;
}

//const PRODART::POSE::atom_type pseudo_N =  PRODART::POSE::atom_type("N'");
//const PRODART::POSE::atom_type pseudo_O =  PRODART::POSE::atom_type("O'");


class ca_pose_meta;

typedef boost::shared_ptr<ca_pose_meta> ca_pose_meta_shared_ptr;
typedef boost::weak_ptr<ca_pose_meta> ca_pose_meta_weak_ptr;

typedef std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> const_atom_shared_ptr_const_atom_shared_ptr_map;

//! class for storing CA pose meta data (pair lists etc) for use during simulations
class ca_pose_meta : public pose_meta_interface {

	friend ca_pose_meta_shared_ptr new_ca_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach,
			const bool auto_load_restraints);

private:
	static bool std_atoms_initialised;
	static PRODART::POSE::atom_type pseudo_N_;
	static PRODART::POSE::atom_type pseudo_O_;

public:

	//these may not be strictly thread safe but sort out later
	static const PRODART::POSE::atom_type pseudo_N(){
		return std_atoms_initialised ? pseudo_N_ : PRODART::POSE::atom_type("N'");
	}
	static const PRODART::POSE::atom_type pseudo_O(){
		return std_atoms_initialised ? pseudo_O_ :  PRODART::POSE::atom_type("O'");
	}


protected:

	void update_pair_lists();

	void update_ca_pair_list();
	void update_no_pair_list();

	ca_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach);


	void initialise_cells();
	//bool assign_atoms_to_cells();

	void initialise_ca_sim_cell();
	bool assign_atoms_to_ca_sim_cell();

	void initialise_no_sim_cell();
	bool assign_atoms_to_no_sim_cell();

	void add_sep_info();
	void add_sep_info(PRODART::POSE::META::nb_ele_vector& pair_list);

	//! TODO to be done by alphabet classifier
	bool quick_add_ca_bonds();


	void set_up_fragments();

	void init();

	void calc_sec_struct();
	void calc_hbond_counts();

	double ca_cutoff;
	double ca_cell_margin;
	PRODART::POSE::META::nb_ele_vector ca_pair_list;
	sim_cell ca_sim_cell;

	PRODART::POSE::atom_shared_ptr_vector N_atoms;
	PRODART::POSE::atom_shared_ptr_vector O_atoms;
	//PRODART::POSE::const_atom_shared_ptr_vector C_i_atoms;
	//PRODART::POSE::const_atom_shared_ptr_vector C_i_p1_atoms;
	const_atom_shared_ptr_ca_hbond_collection_shared_ptr_map ca_hbond_collection_map;
	double NO_cutoff;
	double NO_cell_margin;
	PRODART::POSE::META::nb_ele_vector NO_pair_list;
	sim_cell NO_sim_cell;

	boost::shared_ptr<FRAG_CLASS::fragment_classifier_interface> frag_classifier;

	ca_pose_meta_weak_ptr _this_ca_pose_meta;

	frag4_vector fragments;
	int_vector frag_type_num_vec;



	//
	double pseudo_hb_dist_cutoff;// = 2.6;




private:

public:


	void recalc_pair_lists_dists();


	//! sets CA non-bonded cutoff - this is initialised through prodart_env option key "sim:ca:ca_nb_cutoff" in the constructor and this is the preferred method to change this
	void set_ca_cutoff(const double val);


	PRODART::POSE::META::nb_ele_vector& get_ca_pair_list();
	PRODART::POSE::META::nb_ele_vector& get_NO_pair_list();
	frag4_vector& get_fragments();

	int get_type_frag_num(const int res_num) const{
		return frag_type_num_vec[res_num];
	}

	void set_frag_type_num(const int res_num, const int frag_num) {
		frag_type_num_vec[res_num] = frag_num;
	}

	ca_hbond_collection_shared_ptr get_hbond_collection(const_atom_shared_ptr atm) const{
		return ca_hbond_collection_map.find(atm) != ca_hbond_collection_map.end() ? ca_hbond_collection_map.find(atm)->second : ca_hbond_collection_shared_ptr(); //ca_hbond_collection_map[atm];
	}

	void inactivate_pseudo_NO();
	static void inactivate_pseudo_NO(POSE::pose_shared_ptr this_pose);
	void activate_pseudo_NO();

	poseMetaType get_pose_meta_type() const{
			return pm_ca_pose_meta;
	}


};

ca_pose_meta_shared_ptr new_ca_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach,
		const bool auto_load_restraints = true);


}
}
}

#endif /* CA_POSE_META_H_ */
