//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * ca_pose_meta.cpp
 *
 *  Created on: 5 Mar 2010
 *      Author: jmacdona
 */
#include "ca_pose_meta.h"
#include "frag_classify/fragment_classifier_interface.h"
#include "frag_classify/ca_fragment_classifier.h"

using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;

using namespace std;


namespace PRODART {
namespace POSE {
namespace META {




ca_pose_meta_shared_ptr new_ca_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach,
		const bool auto_load_restraints){
	ca_pose_meta_shared_ptr nmd(new ca_pose_meta(pose_to_attach));
	nmd->_this_ca_pose_meta = nmd;
	nmd->init();
	if (auto_load_restraints){
		POTENTIALS::restraints_store::Instance()->apply_restraints(nmd);
	}
	else {
		cout << "new_ca_pose_meta: INFO: not auto-loading restraints from restraints_store" << endl;
	}
	//npose->_this = nms;
	return nmd;
}


//const PRODART::POSE::atom_type ca_pose_meta::pseudo_N =  PRODART::POSE::atom_type("N'");
//const PRODART::POSE::atom_type ca_pose_meta::pseudo_O =  PRODART::POSE::atom_type("O'");

bool ca_pose_meta::std_atoms_initialised = false;
PRODART::POSE::atom_type ca_pose_meta::pseudo_N_ = PRODART::POSE::atom_type("N'", true);
PRODART::POSE::atom_type ca_pose_meta::pseudo_O_ = PRODART::POSE::atom_type("O'", true);

boost::mutex std_ca_atm_init_mutex;
ca_pose_meta::ca_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach) : pose_meta_interface(pose_to_attach){
	if (!std_atoms_initialised){
		boost::mutex::scoped_lock lock(std_ca_atm_init_mutex);
		pseudo_N_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("N");
		pseudo_O_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("O");
	}
	std_atoms_initialised = true;
}

void ca_pose_meta::init(){



	pseudo_hb_dist_cutoff = PRODART::ENV::get_option_value<double>("sim:ca:pseudo_hb_dist_cutoff"); //4.5;


	this->frag_classifier = FRAG_CLASS::ca_fragment_classifier::MainInstance();
	fragments.resize(0);
	fragments.reserve(pose_->get_residue_count());
	frag_type_num_vec.resize(0);
	frag_type_num_vec.resize(pose_->get_residue_count());
	this->set_up_fragments();
	this->frag_classifier->classify_fragments(pose_, this->_this_ca_pose_meta.lock());

	ca_cell_margin = PRODART::ENV::get_option_value<double>("sim:ca:ca_cell_margin");	//100;
	ca_cutoff = PRODART::ENV::get_option_value<double>("sim:ca:ca_nb_cutoff"); //10;
	ca_pair_list.resize(0);
	ca_pair_list.reserve(pose_->get_residue_count() * 400);
	this->ca_sim_cell.set_cutoff(ca_cutoff);

	this->NO_cell_margin = PRODART::ENV::get_option_value<double>("sim:ca:no_cell_margin");//50;
	this->NO_cutoff = PRODART::ENV::get_option_value<double>("sim:ca:no_nb_cutoff");//10;
	this->NO_pair_list.resize(0);
	this->NO_pair_list.reserve(pose_->get_residue_count() * 400);
	this->NO_sim_cell.set_cutoff(NO_cutoff);

	this->initialise_cells();
	this->update_pair_lists();

	this->calc_hbond_counts();
	this->calc_sec_struct();

	// need this as N' and O' atoms have been added
	this->gradient.resize(0);
	this->gradient.resize(pose_->get_all_atom_count(), vector3d(0,0,0));

	this->updateBondedPairList(this->bond_list);
	this->updateDihedralList(dih_restraint_list);
	this->updateCoordRestraintList(coord_restraint_list);
	this->updateAngleList(this->ang_restraint_list);

	move_margin = PRODART::ENV::get_option_value<double>("sim:ca:move_margin");
	//this->quick_add_ca_bonds();
	make_update_state_store();
	make_has_moved_coord_store();

}


PRODART::POSE::META::nb_ele_vector& ca_pose_meta::get_ca_pair_list(){
	return this->ca_pair_list;
}

PRODART::POSE::META::nb_ele_vector& ca_pose_meta::get_NO_pair_list(){
	return this->NO_pair_list;
}

frag4_vector& ca_pose_meta::get_fragments(){
	return this->fragments;
}

bool ca_pose_meta::assign_atoms_to_ca_sim_cell(){
	ca_sim_cell.clear_cells();
	bool result = true;
	const int resCount = pose_->get_residue_count();
	//const int CA_pos = PRODART::PROT::getRelAtomPosition(PRODART::PROT::CA); //(static_cast<int>(PRODART::PROT::CA)  - (static_cast<int>(PRODART::PROT::UNK_ATOM) +1));

	for (int i = 0; i < resCount; i++){

		const_atom_shared_ptr atm = pose_->get_bb_atom(CA, i);
		if (atm->isActiveAndSet()){
			const bool this_result = ca_sim_cell.add_to_cell(atm);
			result = result && this_result;
			if (result == false) {
				PRINT_EXPR(atm->get_residue()->get_internal_residue_index());
				PRINT_EXPR(atm->get_residue()->get_pdb_residue_index());
				PRINT_EXPR(atm->isActiveAndSet());
				PRINT_EXPR(atm->get_coords());
			}
		}
		if (result == false) return false;
	}

	return result;
}


bool ca_pose_meta::assign_atoms_to_no_sim_cell(){
	NO_sim_cell.clear_cells();
	bool result = true;
	//const int resCount = pose_->get_residue_count();
	//const int CA_pos = PRODART::PROT::getRelAtomPosition(PRODART::PROT::CA); //(static_cast<int>(PRODART::PROT::CA)  - (static_cast<int>(PRODART::PROT::UNK_ATOM) +1));

	PRODART::POSE::atom_shared_ptr_vector::iterator iter;

	for (iter = this->N_atoms.begin(); iter != this->N_atoms.end(); iter++ ){
		//const_atom_shared_ptr atm = pose_->get_bb_atom(CA, i);
		if ((*iter)->isActive()){
			result = result && NO_sim_cell.add_to_cell(*iter);
		}
		if (result == false) return false;
	}
	for (iter = this->O_atoms.begin(); iter != this->O_atoms.end(); iter++ ){
		//const_atom_shared_ptr atm = pose_->get_bb_atom(CA, i);
		if ((*iter)->isActive()){
			result = result && NO_sim_cell.add_to_cell(*iter);
		}
		if (result == false) return false;
	}

	return result;
}

void ca_pose_meta::recalc_pair_lists_dists(){
	//this->updateCA_nb_pair_list_cells(protein);
	this->frag_classifier->classify_fragments(pose_, this->_this_ca_pose_meta.lock());
	update_has_moved_flags();

	if (pair_lists_ok() == false){
		this->refresh_pair_lists = true;
		//cout << "updating...\n";
	}

	if (this->refresh_pair_lists){
		this->update_pair_lists();
	}
	else {
		this->ca_sim_cell.recalcPairListDists(this->ca_pair_list);
		this->NO_sim_cell.recalcPairListDists(this->NO_pair_list);
	}
	refresh_pair_lists = false;
	this->updateBondedPairList(this->bond_list);
	this->updateDihedralList(dih_restraint_list);
	this->updateCoordRestraintList(coord_restraint_list);
	this->updateAngleList(this->ang_restraint_list);
	this->calc_hbond_counts();
	this->calc_sec_struct();
	reset_has_moved_coord_store();
}


void ca_pose_meta::update_ca_pair_list(){
	//ca_pair_list.clear();
	ca_pair_list.resize(0);
	//ca_pair_list.reserve(_pose->get_residue_count() * 400);
	const bool result = this->assign_atoms_to_ca_sim_cell();
	if (result == false){
		cout << "INFO: CAs outside simulation box - reboxing" << endl;
		vector3d highest, lowest;
		PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);
		cout << "lowest:\t" << lowest << endl;
		cout << "highest:\t" << highest << endl;
		this->initialise_ca_sim_cell();
		this->assign_atoms_to_ca_sim_cell();
	}
	//cout << "CA\n" << endl;
	ca_sim_cell.updatePairList_trim(this->ca_pair_list, pose_);
	//this->add_sep_info();
}

void ca_pose_meta::update_no_pair_list(){
	//NO_pair_list.clear();
	NO_pair_list.resize(0);
	//NO_pair_list.reserve(_pose->get_residue_count() * 400);
	const bool result = this->assign_atoms_to_no_sim_cell();
	if (result == false){
		cout << "INFO: NOs outside simulation box - reboxing" << endl;
		vector3d highest, lowest;
		PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);
		cout << "lowest:\t" << lowest << endl;
		cout << "highest:\t" << highest << endl;
		this->initialise_no_sim_cell();
		this->assign_atoms_to_no_sim_cell();
	}
	//cout << "NO\n" << endl;
	NO_sim_cell.updatePairList_trim(this->NO_pair_list, pose_);
}

void ca_pose_meta::update_pair_lists(){
	this->update_ca_pair_list();
	this->update_no_pair_list();
	this->add_sep_info();
	reset_update_state_store();
}

void ca_pose_meta::add_sep_info(){


	this->add_sep_info(this->ca_pair_list);
	this->add_sep_info(this->NO_pair_list);

	/*
	nb_ele_vector::iterator iter;
	for (iter = this->ca_pair_list.begin(); iter != ca_pair_list.end(); iter++){
		if (iter->atom1_ptr->get_chain() == iter->atom2_ptr->get_chain()){
			iter->seq_sep = iter->atom1_ptr->get_residue()->get_internal_residue_index() - iter->atom2_ptr->get_residue()->get_internal_residue_index();
			iter->bond_sep = iter->seq_sep;
		}
		else {
			iter->seq_sep = 99999;
		}
	}
	*/


}

void ca_pose_meta::add_sep_info(PRODART::POSE::META::nb_ele_vector& pair_list){
	nb_ele_vector::iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		iter->res_num1 = iter->atom1_ptr->get_residue()->get_internal_residue_index();
		iter->res_num2 = iter->atom2_ptr->get_residue()->get_internal_residue_index();
		if (iter->atom1_ptr->get_chain() == iter->atom2_ptr->get_chain()){
			iter->seq_sep = std::abs(iter->res_num1 - iter->res_num2);
			iter->bond_sep = iter->seq_sep;
		}
		else {
			iter->seq_sep = 99999;
			iter->bond_sep = iter->seq_sep;
		}
	}
}

void ca_pose_meta::initialise_ca_sim_cell(){
	this->ca_sim_cell.set_cutoff(ca_cutoff);
	vector3d highest, lowest;


	PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);

	//cout << "initialising CA cells\n";

	lowest[0] = lowest[0] - ca_cell_margin;
	lowest[1] = lowest[1] - ca_cell_margin;
	lowest[2] = lowest[2] - ca_cell_margin;

	PRINT_EXPR(lowest);

	//cout << "lowest:\t" << lowest << endl;

	highest[0] = highest[0] + ca_cell_margin;
	highest[1] = highest[1] + ca_cell_margin;
	highest[2] = highest[2] + ca_cell_margin;

	PRINT_EXPR(highest);

	//cout << "highest:\t" << highest << endl;

	this->ca_sim_cell.set_x_limits(lowest[0], highest[0]);
	this->ca_sim_cell.set_y_limits(lowest[1], highest[1]);
	this->ca_sim_cell.set_z_limits(lowest[2], highest[2]);

	this->ca_sim_cell.set_num_atoms(pose_->get_residue_count());

	//this->ca_sim_cell.set_cutoff(ca_cutoff);

	this->ca_sim_cell.create_cells();
	this->ca_sim_cell.create_cell_neighbour_vec();

}

void ca_pose_meta::initialise_no_sim_cell(){
	this->NO_sim_cell.set_cutoff(NO_cutoff);
	vector3d highest, lowest;


	PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);

	//cout << "initialising NO cells\n";

	lowest[0] = lowest[0] - NO_cell_margin;
	lowest[1] = lowest[1] - NO_cell_margin;
	lowest[2] = lowest[2] - NO_cell_margin;

	//cout << "lowest:\t" << lowest << endl;

	highest[0] = highest[0] + NO_cell_margin;
	highest[1] = highest[1] + NO_cell_margin;
	highest[2] = highest[2] + NO_cell_margin;

	//cout << "highest:\t" << highest << endl;

	this->NO_sim_cell.set_x_limits(lowest[0], highest[0]);
	this->NO_sim_cell.set_y_limits(lowest[1], highest[1]);
	this->NO_sim_cell.set_z_limits(lowest[2], highest[2]);

	this->NO_sim_cell.set_num_atoms(pose_->get_residue_count()*2);

	//this->NO_sim_cell.set_cutoff(NO_cutoff);

	this->NO_sim_cell.create_cells();
	this->NO_sim_cell.create_cell_neighbour_vec();
}

class ca_pose_meta_defunct_option_exception: public std::exception{
  virtual const char* what() const throw()
  {
    return "ERROR: set this through the prodart_env option manager instead";
  }
};

void ca_pose_meta::set_ca_cutoff(const double val){
	this->ca_cutoff = val;
}


bool ca_pose_meta::quick_add_ca_bonds(){

	const int resCount = this->pose_->get_residue_count();
	// TODO should this be cleared???
	this->bond_list.resize(0);
	this->bond_list.reserve(resCount + 10);

	if (resCount <= 0){
		return false;
	}

	residue_shared_ptr prev_res;

	for (int i = 0; i < resCount-1; i++){
		residue_shared_ptr this_res = pose_->get_residue(i);
		if (prev_res && this_res){

			if (prev_res->get_chain() == this_res->get_chain()){

				atom_shared_ptr at1 = pose_->get_bb_atom(CA, i);
				atom_shared_ptr at2 = pose_->get_bb_atom(CA, i+1);
				if (at1 && at2){
					// add bond
					bonded_pair_element ele;
					ele.atom1_ptr = at1;
					ele.atom2_ptr = at2;
					ele.atom1 = at1->get_seq_num();
					ele.atom2 = at2->get_seq_num();
					ele.equilibrium_dist = 3.8;
					ele.half_bond_const = 10.0;
					this->bond_list.push_back(ele);
				}


			}

		}
		prev_res = this_res;
	}



	return true;
}


void ca_pose_meta::initialise_cells(){
	this->initialise_ca_sim_cell();
	this->initialise_no_sim_cell();
}

/*
bool ca_pose_meta::assign_atoms_to_cells(){
	return this->assign_atoms_to_ca_sim_cell();
}
*/



// Question: Should I be checking for Active or Set and doesn't only one of 4 need to be Active?
// Answer: in fact I am doing the right thing here so no need for change
void ca_pose_meta::set_up_fragments(){

	this->fragments.resize(0);
	this->fragments.reserve(pose_->get_residue_count());

	N_atoms.resize(0);
	N_atoms.reserve(pose_->get_residue_count());

	O_atoms.resize(0);
	O_atoms.reserve(pose_->get_residue_count());

	ca_hbond_collection_map.clear();

	/*
	this->C_i_atoms.clear();
	C_i_atoms.reserve(pose_->get_residue_count());

	this->C_i_p1_atoms.clear();
	C_i_p1_atoms.reserve(pose_->get_residue_count());
	*/

	if ( pose_->get_residue_count() < 4){
		return;
	}

	//PRODART::POSE::const_atom_shared_ptr atoms[4];

	class frag4_element ele;
	//ele.CA_atoms[0] = _pose->get_bb_atom(CA, 0);
	ele.CA_atoms[1] = pose_->get_bb_atom(CA, 0);
	ele.CA_atoms[2] = pose_->get_bb_atom(CA, 1);
	ele.CA_atoms[3] = pose_->get_bb_atom(CA, 2);

	if (ele.CA_atoms[1]->get_chain() == ele.CA_atoms[2]->get_chain()
			&& ele.CA_atoms[2]->get_chain() == ele.CA_atoms[3]->get_chain()
			&& ele.CA_atoms[1]->isActive()
			&& ele.CA_atoms[2]->isActive()
			&& ele.CA_atoms[3]->isActive()){
		PRODART::POSE::atom_shared_ptr N_prime = pose_->add_new_atom(vector3d(0,0,0), pseudo_N(), 0);
		PRODART::POSE::atom_shared_ptr O_prime = pose_->add_new_atom(vector3d(0,0,0), pseudo_O(), 0);
		N_atoms.push_back(N_prime);
		O_atoms.push_back(O_prime);
		ca_hbond_collection_shared_ptr hb_ele(new ca_hbond_collection());
		hb_ele->N_i = N_prime;
		hb_ele->O_i = O_prime;
		hb_ele->CA_i = ele.CA_atoms[1];
		hb_ele->CA_i_p1 = ele.CA_atoms[2];
		this->ca_hbond_collection_map[N_prime] = hb_ele;
		this->ca_hbond_collection_map[O_prime] = hb_ele;
		this->ca_hbond_collection_map[hb_ele->CA_i] = hb_ele;
		//this->ca_hbond_collection_map[hb_ele->CA_i_p1] = hb_ele;
	}

	for (int i = 3; i < pose_->get_residue_count(); i++){
		ele.CA_atoms[0] = 	ele.CA_atoms[1];
		ele.CA_atoms[1] = 	ele.CA_atoms[2];
		ele.CA_atoms[2] = 	ele.CA_atoms[3];
		ele.CA_atoms[3] = pose_->get_bb_atom(CA, i);

		if (ele.CA_atoms[0]->isActive()
				&& ele.CA_atoms[1]->isActive()
				&& ele.CA_atoms[2]->isActive()
				&& ele.CA_atoms[3]->isActive()){

			if (ele.CA_atoms[0]->get_chain() == ele.CA_atoms[1]->get_chain()
					&& ele.CA_atoms[1]->get_chain() == ele.CA_atoms[2]->get_chain()){
				PRODART::POSE::atom_shared_ptr N_prime = pose_->add_new_atom(vector3d(0,0,0), pseudo_N(), i-2);
				PRODART::POSE::atom_shared_ptr O_prime = pose_->add_new_atom(vector3d(0,0,0), pseudo_O(), i-2);
				N_atoms.push_back(N_prime);
				O_atoms.push_back(O_prime);
				ca_hbond_collection_shared_ptr hb_ele(new ca_hbond_collection());
				hb_ele->N_i = N_prime;
				hb_ele->O_i = O_prime;
				hb_ele->CA_i = ele.CA_atoms[1];
				hb_ele->CA_i_p1 = ele.CA_atoms[2];
				this->ca_hbond_collection_map[N_prime] = hb_ele;
				this->ca_hbond_collection_map[O_prime] = hb_ele;
				this->ca_hbond_collection_map[hb_ele->CA_i] = hb_ele;
				//this->ca_hbond_collection_map[hb_ele->CA_i_p1] = hb_ele;
			}

			if (ele.CA_atoms[0]->get_chain() == ele.CA_atoms[1]->get_chain()
					&& ele.CA_atoms[1]->get_chain() == ele.CA_atoms[2]->get_chain()
					&& ele.CA_atoms[2]->get_chain() == ele.CA_atoms[3]->get_chain()){
				if (ele.CA_atoms[0]->get_residue()->get_prev_residue()
						&& ele.CA_atoms[3]->get_residue()->get_next_residue() ){
					ele.frag_pos = fragMIDDLE;
				}
				else if ( !ele.CA_atoms[0]->get_residue()->get_prev_residue()
						&& !ele.CA_atoms[3]->get_residue()->get_next_residue()){
					ele.frag_pos = fragNandCTERM;
				}
				else if (!ele.CA_atoms[0]->get_residue()->get_prev_residue()){
					ele.frag_pos = fragNTERM;
				}
				else if (!ele.CA_atoms[3]->get_residue()->get_next_residue()) {
					ele.frag_pos = fragCTERM;
				}
				else {
					//some error!
					std::cerr << "\nWTF: SOME ERROR SETTING UP FRAGMENTS!!!!!!!!!!\n" << std::endl;
					ele.frag_pos = fragMIDDLE;
				}

				ele.frag_type_num = 0;
				ele.frag_ss_class = PRODART::POSE::ss4_OTHER;
				ele.CA_1_residue_number = i - 2;
				this->fragments.push_back(ele);
			}

		}
	}




}


void ca_pose_meta::calc_hbond_counts(){

	strand_hbond_count.resize(0);
	strand_hbond_count.resize(pose_->get_residue_count(), 0);
	helix_hbond_count.resize(0);
	helix_hbond_count.resize(pose_->get_residue_count(), 0);

	PRODART::POSE::META::nb_ele_vector& pair_list = this->get_NO_pair_list();
	PRODART::POSE::four_state_sec_struct_vector& conf_class = this->get_conf_class();
	PRODART::POSE::META::nb_ele_vector::iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		//cout << "\ndb: 1 ";
		if (iter->dist < this->pseudo_hb_dist_cutoff){
			//cout << "db: 2 ";
			const int seq_sep = iter->seq_sep;
			PRODART::POSE::four_state_sec_struct ss1 = conf_class[iter->res_num1];
			PRODART::POSE::four_state_sec_struct ss2 = conf_class[iter->res_num2];
			//cout << "db: " << iter->atype1.get_label() << " " << iter->atype2.get_label() << " ";
			if ((iter->atype1 == pseudo_N() &&  iter->atype2 == pseudo_O())
					|| (iter->atype1 == pseudo_O() &&  iter->atype2 == pseudo_N()) ){


				if (seq_sep == 3
						&& (ss1 == PRODART::POSE::ss4_HELIX // secs_HELIX
						|| ss2 == PRODART::POSE::ss4_HELIX)){
					helix_hbond_count[iter->res_num1]++;
					helix_hbond_count[iter->res_num2]++;
				}
				else if (seq_sep >= 3
						&& (ss1 == PRODART::POSE::ss4_STRAND
						|| ss2 == PRODART::POSE::ss4_STRAND)){
					strand_hbond_count[iter->res_num1]++;
					strand_hbond_count[iter->res_num2]++;
				}

			}
		}
	}
}

void ca_pose_meta::calc_sec_struct(){
	const int resCount = pose_->get_residue_count();
	this->residue_sec_struct.resize(0);
	this->residue_sec_struct.resize(resCount, PRODART::POSE::ss3_OTHER);

	for (int i = 1; i < resCount-1; i++){
		if (helix_hbond_count[i] > 0
				&& residue_conformation_class[i] == PRODART::POSE::ss4_HELIX){
			residue_sec_struct[i] = PRODART::POSE::ss3_HELIX;//secs_HELIX;
		}
		else if (strand_hbond_count[i] > 0
				&& residue_conformation_class[i] == PRODART::POSE::ss4_STRAND){
			residue_sec_struct[i] = PRODART::POSE::ss3_STRAND;// secs_STRAND;
		}
		else if ((residue_conformation_class[i] == PRODART::POSE::ss4_STRAND
				)
				&& strand_hbond_count[i-1] > 0
				&& strand_hbond_count[i+1] > 0){
			residue_sec_struct[i] = PRODART::POSE::ss3_STRAND;// secs_STRAND;
		}
		else if ((residue_conformation_class[i] == PRODART::POSE::ss4_HELIX
				|| helix_hbond_count[i] > 0)
				&& helix_hbond_count[i-1] > 0
				&& helix_hbond_count[i+1] > 0){
			residue_sec_struct[i] = PRODART::POSE::ss3_HELIX;// secs_HELIX;
		}
		else {
			residue_sec_struct[i] = PRODART::POSE::ss3_OTHER;// secs_OTHER;
		}
	}

	PRODART::POSE::three_state_sec_struct prev_val = PRODART::POSE::ss3_OTHER;//  secs_OTHER;
	PRODART::POSE::three_state_sec_struct next_val = PRODART::POSE::ss3_OTHER;//  secs_OTHER;
	for (int i = 0; i < resCount - 1; i++){
		next_val = residue_sec_struct[i+1];
		if (residue_sec_struct[i] != PRODART::POSE::ss3_OTHER){
			if (residue_sec_struct[i] != prev_val
			      && residue_sec_struct[i] != next_val){
				residue_sec_struct[i] = PRODART::POSE::ss3_OTHER;
			}
		}
		prev_val = residue_sec_struct[i];
	}
}


void ca_pose_meta::inactivate_pseudo_NO(){

	const int resCount = pose_->get_residue_count();


	for (int i = 0; i < resCount; i++){

		atom_shared_ptr psN = pose_->get_atom(pseudo_N(), i);
		atom_shared_ptr psO = pose_->get_atom(pseudo_O(), i);

		if (psN){
			psN->setActive(false);
		}
		if (psO){
			psO->setActive(false);
		}

	}

}

void ca_pose_meta::inactivate_pseudo_NO(POSE::pose_shared_ptr this_pose){
	const int resCount = this_pose->get_residue_count();


	for (int i = 0; i < resCount; i++){

		atom_shared_ptr psN = this_pose->get_atom(pseudo_N(), i);
		atom_shared_ptr psO = this_pose->get_atom(pseudo_O(), i);

		if (psN){
			psN->setActive(false);
		}
		if (psO){
			psO->setActive(false);
		}

	}
}

void ca_pose_meta::activate_pseudo_NO(){

	const int resCount = pose_->get_residue_count();


	for (int i = 0; i < resCount; i++){

		atom_shared_ptr psN = pose_->get_atom(pseudo_N(), i);
		atom_shared_ptr psO = pose_->get_atom(pseudo_O(), i);

		if (psN){
			psN->setActive(true);
		}
		if (psO){
			psO->setActive(true);
		}

	}

}

}
}
}




