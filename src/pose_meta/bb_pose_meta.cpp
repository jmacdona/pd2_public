/*
 * bb_pose_meta.cpp
 *
 *  Created on: 22 Oct 2010
 *      Author: jmacdona
 */

#include "bb_pose_meta.h"

using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;

using namespace std;


namespace PRODART {
namespace POSE {
namespace META {

bb_pose_meta_shared_ptr new_bb_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach){
	bb_pose_meta_shared_ptr nmd(new bb_pose_meta(pose_to_attach));
	nmd->_this_bb_pose_meta = nmd;
	nmd->init();
	POTENTIALS::restraints_store::Instance()->apply_restraints(nmd);
	//npose->_this = nms;
	return nmd;
}

const double bb_pose_meta::num_bb_atoms = 6;

bool bb_pose_meta::std_atoms_initialised = false;
PRODART::POSE::atom_type bb_pose_meta::bb_aty_N_ =  PRODART::POSE::atom_type("N", true);
PRODART::POSE::atom_type bb_pose_meta::bb_aty_O_ =  PRODART::POSE::atom_type("O", true);
PRODART::POSE::atom_type bb_pose_meta::bb_aty_H_ =  PRODART::POSE::atom_type("H", true);
PRODART::POSE::atom_type bb_pose_meta::bb_aty_CA_ =  PRODART::POSE::atom_type("CA", true);
PRODART::POSE::atom_type bb_pose_meta::bb_aty_CB_ =  PRODART::POSE::atom_type("CB", true);
PRODART::POSE::atom_type bb_pose_meta::bb_aty_C_ =  PRODART::POSE::atom_type("C", true);

//const prot_backbone_map* const bb_pose_meta::bb_map(prot_backbone_map::Instance());

boost::mutex std_atm_init_mutex;
bb_pose_meta::bb_pose_meta(PRODART::POSE::pose_shared_ptr pose_to_attach) : pose_meta_interface(pose_to_attach), bb_map(prot_backbone_map::Instance()){

	if (!std_atoms_initialised){
		boost::mutex::scoped_lock lock(std_atm_init_mutex);
		bb_aty_N_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("N");
		bb_aty_O_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("O");
		bb_aty_H_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("H");
		bb_aty_CA_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("CA");
		bb_aty_CB_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("CB");
		bb_aty_C_.auto_set_ff_atom_id();// =  PRODART::POSE::atom_type("C");
	}
	std_atoms_initialised = true;

}



void bb_pose_meta::recalc_pair_lists_dists(){

	update_has_moved_flags();
	if (pair_lists_ok() == false){
		this->refresh_pair_lists = true;
		//cout << "updating...\n";
	}

	if (this->refresh_pair_lists){
		this->update_pair_lists();
	}
	else {
		this->all_bb_sim_cell.recalcPairListDists(this->all_bb_pair_list);
		this->H_O_sim_cell.recalcPairListDists(this->H_O_pair_list);
	}
	this->updateBondedPairList(this->bb_bond_list);
	this->updateBondedPairList(this->bond_list);
	this->updateDihedralList(dih_restraint_list);
	this->updateCoordRestraintList(coord_restraint_list);
	this->updateAngleList(this->bb_angle_list);
	this->updateAngleList(this->ang_restraint_list);
	this->updateDihedralList(this->bb_dih_angle_list);
	this->update_phi_psi_omega();
	this->updateNonBondedPairList(this->bonded_14_list);
	this->updateNonBondedPairList(this->bonded_15_list);

	this->calc_hbond_counts();
	this->calc_sec_struct();
	reset_has_moved_coord_store();
	refresh_pair_lists = false;
}

PRODART::POSE::META::nb_ele_vector& bb_pose_meta::get_all_bb_pair_list(){
	return this->all_bb_pair_list;
}

PRODART::POSE::META::nb_ele_vector& bb_pose_meta::get_H_O_pair_list(){
	return this->H_O_pair_list;
}

PRODART::POSE::META::nb_ele_vector& bb_pose_meta::get_bonded_14_list(){
	return this->bonded_14_list;
}

PRODART::POSE::META::nb_ele_vector& bb_pose_meta::get_bonded_15_list(){
	return this->bonded_15_list;
}

PRODART::POSE::META::bonded_ele_vector& bb_pose_meta::get_bb_bond_list(){
	return this->bb_bond_list;
}

PRODART::POSE::META::angle_ele_vector& bb_pose_meta::get_bb_angle_list(){
	return this->bb_angle_list;
}

PRODART::POSE::META::dihedral_ele_vector& bb_pose_meta::get_bb_dih_angle_list(){
	return this->bb_dih_angle_list;
}

void bb_pose_meta::init(){

	all_bb_cutoff = PRODART::ENV::get_option_value<double>("sim:bb:all_bb_nb_cutoff");
	all_bb_cell_margin = PRODART::ENV::get_option_value<double>("sim:bb:all_bb_cell_margin");
	H_O_cutoff = PRODART::ENV::get_option_value<double>("sim:bb:ho_nb_cutoff");
	H_O_cell_margin = PRODART::ENV::get_option_value<double>("sim:bb:ho_cell_margin");
	hb_dist_cutoff = PRODART::ENV::get_option_value<double>("sim:bb:hb_dist_cutoff");

	this->regularise_atom_names();

	this->initialise_cells();
	this->update_pair_lists();

	this->createBondedPairList();
	this->createAngleList();
	this->createDihedralList();
	this->createBonded14_15Lists();

	this->updateBondedPairList(this->bond_list);
	this->updateDihedralList(dih_restraint_list);
	this->updateCoordRestraintList(coord_restraint_list);
	this->updateBondedPairList(this->bb_bond_list);
	this->updateAngleList(this->bb_angle_list);
	this->updateAngleList(this->ang_restraint_list);
	this->updateDihedralList(this->bb_dih_angle_list);
	this->update_phi_psi_omega();
	this->updateNonBondedPairList(this->bonded_14_list);
	this->updateNonBondedPairList(this->bonded_15_list);

	this->calc_hbond_counts();
	this->calc_sec_struct();

	move_margin = PRODART::ENV::get_option_value<double>("sim:bb:move_margin");
	make_update_state_store();
	make_has_moved_coord_store();

}

void bb_pose_meta::update_pair_lists(){

	this->update_all_bb_pair_list();
	this->update_H_O_pair_list();
	this->add_sep_info();
	reset_update_state_store();
}

void bb_pose_meta::initialise_cells(){
	this->initialise_all_bb_cell();
	this->initialise_H_O_cell();
}
void bb_pose_meta::initialise_all_bb_cell(){
	this->all_bb_sim_cell.set_cutoff(this->all_bb_cutoff);
	vector3d highest, lowest;


	PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);

	//cout << "initialising all BB cells\n";

	lowest[0] = lowest[0] - all_bb_cell_margin;
	lowest[1] = lowest[1] - all_bb_cell_margin;
	lowest[2] = lowest[2] - all_bb_cell_margin;

	//cout << "lowest:\t" << lowest << endl;

	highest[0] = highest[0] + all_bb_cell_margin;
	highest[1] = highest[1] + all_bb_cell_margin;
	highest[2] = highest[2] + all_bb_cell_margin;

	//cout << "highest:\t" << highest << endl;

	this->all_bb_sim_cell.set_x_limits(lowest[0], highest[0]);
	this->all_bb_sim_cell.set_y_limits(lowest[1], highest[1]);
	this->all_bb_sim_cell.set_z_limits(lowest[2], highest[2]);

	this->all_bb_sim_cell.set_num_atoms(pose_->get_residue_count() * num_bb_atoms);

	//this->all_bb_sim_cell.set_cutoff(ca_cutoff);

	this->all_bb_sim_cell.create_cells();
	this->all_bb_sim_cell.create_cell_neighbour_vec();
}
void bb_pose_meta::initialise_H_O_cell(){
	this->H_O_sim_cell.set_cutoff(H_O_cutoff);
	vector3d highest, lowest;


	PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);

	//cout << "initialising HO cells\n";

	lowest[0] = lowest[0] - H_O_cell_margin;
	lowest[1] = lowest[1] - H_O_cell_margin;
	lowest[2] = lowest[2] - H_O_cell_margin;

	//cout << "lowest:\t" << lowest << endl;

	highest[0] = highest[0] + H_O_cell_margin;
	highest[1] = highest[1] + H_O_cell_margin;
	highest[2] = highest[2] + H_O_cell_margin;

	//cout << "highest:\t" << highest << endl;

	this->H_O_sim_cell.set_x_limits(lowest[0], highest[0]);
	this->H_O_sim_cell.set_y_limits(lowest[1], highest[1]);
	this->H_O_sim_cell.set_z_limits(lowest[2], highest[2]);

	this->H_O_sim_cell.set_num_atoms(pose_->get_residue_count()*2);

	//this->H_O_sim_cell.set_cutoff(H_O_cutoff);

	this->H_O_sim_cell.create_cells();
	this->H_O_sim_cell.create_cell_neighbour_vec();
}


bool bb_pose_meta::assign_atoms_to_all_bb_sim_cell(){
	all_bb_sim_cell.clear_cells();
	bool result = true;
	const int resCount = pose_->get_residue_count();
	//const int CA_pos = PRODART::PROT::getRelAtomPosition(PRODART::PROT::CA); //(static_cast<int>(PRODART::PROT::CA)  - (static_cast<int>(PRODART::PROT::UNK_ATOM) +1));

	for (int i = 0; i < resCount; i++){

		const_atom_shared_ptr atm = pose_->get_bb_atom(N, i);
		if (atm->isActive()){
			result = result && all_bb_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
		atm = pose_->get_bb_atom(CA, i);
		if (atm->isActive()){
			result = result && all_bb_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
		atm = pose_->get_bb_atom(C, i);
		if (atm->isActive()){
			result = result && all_bb_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
		atm = pose_->get_bb_atom(O, i);
		if (atm->isActive()){
			result = result && all_bb_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
		/*
		atm = pose_->get_bb_atom(H, i);
		if (atm->isActive()){
			result = result && all_bb_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
		*/
		atm = pose_->get_bb_atom(CB, i);
		if (atm->isActive()){
			result = result && all_bb_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
	}

	return result;
}

bool bb_pose_meta::assign_atoms_to_H_O_sim_cell(){
	H_O_sim_cell.clear_cells();
	bool result = true;
	const int resCount = pose_->get_residue_count();
	//const int CA_pos = PRODART::PROT::getRelAtomPosition(PRODART::PROT::CA); //(static_cast<int>(PRODART::PROT::CA)  - (static_cast<int>(PRODART::PROT::UNK_ATOM) +1));

	for (int i = 0; i < resCount; i++){
		const_atom_shared_ptr atm = pose_->get_bb_atom(H, i);
		if (atm->isActive()){
			result = result && H_O_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
		atm = pose_->get_bb_atom(O, i);
		if (atm->isActive()){
			result = result && H_O_sim_cell.add_to_cell(atm);
		}
		if (result == false) return false;
	}

	return result;

}

void bb_pose_meta::update_all_bb_pair_list(){
	all_bb_pair_list.resize(0);
	const bool result = this->assign_atoms_to_all_bb_sim_cell();
	if (result == false){
		cout << "INFO: BB atoms outside simulation box - reboxing" << endl;
		vector3d highest, lowest;
		PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);
		cout << "lowest:\t" << lowest << endl;
		cout << "highest:\t" << highest << endl;
		this->initialise_all_bb_cell();
		this->assign_atoms_to_all_bb_sim_cell();
	}
	all_bb_sim_cell.updatePairList_trim(this->all_bb_pair_list, pose_);
}

void bb_pose_meta::update_H_O_pair_list(){
	H_O_pair_list.resize(0);
	//H_O_pair_list.reserve(_pose->get_residue_count() * 400);
	const bool result = this->assign_atoms_to_H_O_sim_cell();
	if (result == false){
		cout << "INFO: HOs outside simulation box - reboxing" << endl;
		vector3d highest, lowest;
		PRODART::POSE_UTILS::get_ca_bounding_box(this->pose_, lowest, highest);
		cout << "lowest:\t" << lowest << endl;
		cout << "highest:\t" << highest << endl;
		this->initialise_H_O_cell();
		this->assign_atoms_to_H_O_sim_cell();
	}
	//cout << "NO\n" << endl;
	H_O_sim_cell.updatePairList_trim(this->H_O_pair_list, pose_);
}



void bb_pose_meta::add_sep_info(){
	this->add_sep_info(this->all_bb_pair_list);
	this->add_sep_info(this->H_O_pair_list);
}
void bb_pose_meta::add_sep_info(PRODART::POSE::META::nb_ele_vector& pair_list){
	nb_ele_vector::iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		iter->res_num1 = iter->atom1_ptr->get_residue()->get_internal_residue_index();
		iter->res_num2 = iter->atom2_ptr->get_residue()->get_internal_residue_index();
		if (iter->atom1_ptr->get_chain() == iter->atom2_ptr->get_chain()){
			iter->seq_sep = std::abs(iter->res_num1 - iter->res_num2);
			//iter->bond_sep = iter->seq_sep;
			if (iter->res_num2 >= iter->res_num1){
				iter->bond_sep = bb_map->get_bond_sep(iter->seq_sep, iter->atype1, iter->atype2);
			}
			else {
				iter->bond_sep = bb_map->get_bond_sep(iter->seq_sep, iter->atype2, iter->atype1);
			}
		}
		else {
			iter->seq_sep = 99999;
			iter->bond_sep = iter->seq_sep;
		}
	}
}


void bb_pose_meta::add_bb_bonded_pair(PRODART::POSE::atom_shared_ptr at1,
		PRODART::POSE::atom_shared_ptr at2,
		const double eq_dist,
		const double half_bond_const){


	if (at1 && at2){
		bonded_pair_element ele;
		ele.atom1_ptr = at1;
		ele.atom2_ptr = at2;
		ele.atom1 = at1->get_seq_num();
		ele.atom2 = at2->get_seq_num();
		ele.equilibrium_dist = eq_dist;
		ele.half_bond_const = half_bond_const;
		this->bb_bond_list.push_back(ele);

#ifndef NDEBUG
		if ((!at1->isActive()) || (!at2->isActive())){
			std::cerr << "WARNING: bb_pose_meta: add_bb_bonded_pair: inactive backbone atoms" << std::endl;
		}
#endif


	}
	else{
#ifndef NDEBUG
		std::cerr << "WARNING: bb_pose_meta: add_bb_bonded_pair: missing backbone atoms" << std::endl;
#endif
	}

}


void bb_pose_meta::createBondedPairList(){

	const double scale_to_joules = 4184; //scale kcal to J
	const double R = 8.314472; //J K-1 mol-1
	const double assumed_temperature = 200.0;
	const double scale_factor = (scale_to_joules / (R * assumed_temperature)) / 2.0;
	//const double e_r = 5.0;		//80.0 in water
	//cout << "scale_factor: " << scale_factor << endl;

	this->bb_bond_list.resize(0);
	this->bb_bond_list.reserve(pose_->get_residue_count() * 50);

	POSE::residue_shared_ptr res_0, res_1;
	POSE::chain_shared_ptr chain_0, chain_1;
	//cout << "1\n";
	const int resCount = pose_->get_residue_count();
	res_1 = pose_->get_residue(0);
	chain_1 = res_1->get_chain();
	for (int i = 0; i < resCount ; i++){
		//cout << "2";
		res_0 = res_1;
		chain_0 = chain_1;
		if (i + 1 < resCount ) {
			res_1 = pose_->get_residue(i+1);	//protein.getResidue( i + 1 );
		}
		else {
			res_1 = residue_shared_ptr();
		}

		if (res_1){
			chain_1 = res_1->get_chain();
		}
		else {
			chain_1 = chain_shared_ptr();
		}
		atom_shared_ptr atI_H = pose_->get_bb_atom(H, i);//protein.getAtomIndex( i, PRODART::PROT::H );
		atom_shared_ptr atI_N = pose_->get_bb_atom( N, i );//protein.getAtomIndex( i, PRODART::PROT::N );
		atom_shared_ptr atI_CA = pose_->get_bb_atom(CA, i);//protein.getAtomIndex( i, PRODART::PROT::CA );
		atom_shared_ptr atI_CB = pose_->get_bb_atom(CB, i);//protein.getAtomIndex( i, PRODART::PROT::CB );
		atom_shared_ptr atI_C = pose_->get_bb_atom(C, i);//protein.getAtomIndex( i, PRODART::PROT::C );
		atom_shared_ptr atI_O = pose_->get_bb_atom(O, i);//protein.getAtomIndex( i, PRODART::PROT::O );


		/*
		bondLengthParams[CA_C].k = 317.0;
		bondLengthParams[CA_C].l_0 = 1.5220;

		bondLengthParams[C_O].k = 570.0;
		bondLengthParams[C_O].l_0 = 1.2290;

		bondLengthParams[C_N].k = 490.0;
		bondLengthParams[C_N].l_0 = 1.3350;

		bondLengthParams[N_CA].k = 337.0;
		bondLengthParams[N_CA].l_0 = 1.4490;

		bondLengthParams[CA_CB].k = 260.0;
		bondLengthParams[CA_CB].l_0 = 1.5260;

		bondLengthParams[N_H].k = 434.0;
		bondLengthParams[N_H].l_0 = 1.0100;
		*/

		// TODO remove hard coding
		this->add_bb_bonded_pair(atI_H, atI_N, 1.0100, scale_factor * 434.0);
		this->add_bb_bonded_pair(atI_CA, atI_N, 1.4490, scale_factor * 337.0);
		this->add_bb_bonded_pair(atI_CA, atI_CB, 1.5260, scale_factor * 260.0);
		this->add_bb_bonded_pair(atI_CA, atI_C, 1.5220, scale_factor * 317.0);
		this->add_bb_bonded_pair(atI_C, atI_O, 1.2290, scale_factor * 570.0 );

		if (chain_0 == chain_1 && res_1 && chain_1){
			atom_shared_ptr atI_N_p1 = pose_->get_bb_atom(N, i+1);//protein.getAtomIndex( i+1, PRODART::PROT::N );
			this->add_bb_bonded_pair(atI_C, atI_N_p1, 1.3350,  scale_factor * 490.0);
		}


	}
}






void bb_pose_meta::add_bb_angle(PRODART::POSE::atom_shared_ptr at1,
		PRODART::POSE::atom_shared_ptr at2,
		PRODART::POSE::atom_shared_ptr at3,
		const double eq_angle,
		const double half_angle_const){
	if (at1 && at2 && at3){
		angle_element ele;
		ele.atom1_ptr = at1;
		ele.atom2_ptr = at2;
		ele.atom3_ptr = at3;
		ele.atom1 = at1->get_seq_num();
		ele.atom2 = at2->get_seq_num();
		ele.atom3 = at3->get_seq_num();
		ele.res_num1 = at1->get_residue()->get_internal_residue_index();
		ele.res_num2 = at2->get_residue()->get_internal_residue_index();
		ele.res_num3 = at3->get_residue()->get_internal_residue_index();
		ele.equilibrium_angle = eq_angle;
		ele.half_angle_const = half_angle_const;
		this->bb_angle_list.push_back(ele);


#ifndef NDEBUG
		if ((!at1->isActive()) || (!at2->isActive())|| (!at3->isActive())){
			std::cerr << "WARNING: bb_pose_meta: add_bb_angle: inactive backbone atoms" << std::endl;
		}
#endif

	}
	else{
		//std::cerr << "WARNING: bb_pose_meta: add_bb_angle: missing backbone atoms" << std::endl;
	}
}



void bb_pose_meta::createAngleList(){

	const double scale_to_joules = 4184; //scale kcal to J
	const double R = 8.314472; //J K-1 mol-1
	const double assumed_temperature = 200.0;
	const double scale_factor = (scale_to_joules / (R * assumed_temperature)) / 2.0;
	//const double e_r = 5.0;		//80.0 in water
	//cout << "scale_factor: " << scale_factor << endl;


	bb_angle_list.resize(0);
	bb_angle_list.reserve(pose_->get_residue_count() * 50);
	POSE::residue_shared_ptr res_0, res_p1, res_m1;
	POSE::chain_shared_ptr chain_0, chain_1, chain_m1;
	//cout << "1\n";

	const int resCount = pose_->get_residue_count();
	res_m1 = residue_shared_ptr();
	res_p1 = pose_->get_residue(0);
	chain_1 = res_p1->get_chain();
	for (int i = 0; i < resCount ; i++){
		//cout << "2";
		res_m1 = res_0;
		res_0 = res_p1;
		chain_m1 = chain_0;
		chain_0 = chain_1;
		if (i+1 < resCount){
			res_p1 = pose_->get_residue( i + 1 );
		}
		else {
			res_p1 = residue_shared_ptr();
		}

		if (res_p1 ){
			chain_1 = res_p1->get_chain();
		}
		else {
			chain_1 = chain_shared_ptr();
		}

		POSE::atom_shared_ptr ca_0, c_0, n_0, h_0, cb_0, o_0, c_m1, n_p1;

		if (i-1 > 0){
			c_m1 = pose_->get_bb_atom(POSE::C, i-1 );
		}
		else {
			c_m1 = atom_shared_ptr();
		}
		h_0 = pose_->get_bb_atom( POSE::H, i );
		n_0 = pose_->get_bb_atom( POSE::N, i);
		ca_0 = pose_->get_bb_atom(POSE::CA, i);
		cb_0 = pose_->get_bb_atom(POSE::CB, i);
		c_0 = pose_->get_bb_atom(POSE::C, i);
		o_0 = pose_->get_bb_atom(POSE::O ,i);
		if (i+1 < resCount){
			n_p1 = pose_->get_bb_atom( POSE::N, i+1 );
		}
		else {
			n_p1 = atom_shared_ptr();
		}



		/*
		bondAngleParams[N_C_O].k = 80.0;
		bondAngleParams[N_C_O].l_0 = 122.9*(PI / 180.0);

		bondAngleParams[CA_C_O].k = 80.0;
		bondAngleParams[CA_C_O].l_0 = 120.4*(PI / 180.0);

		bondAngleParams[CA_C_N].k = 70.0;
		bondAngleParams[CA_C_N].l_0 = 116.6*(PI / 180.0);

		bondAngleParams[C_N_H].k = 35.0;
		bondAngleParams[C_N_H].l_0 = 119.8*(PI / 180.0);

		bondAngleParams[C_N_CA].k = 50.0;
		bondAngleParams[C_N_CA].l_0 = 121.9*(PI / 180.0);

		bondAngleParams[CA_N_H].k = 38.0;
		bondAngleParams[CA_N_H].l_0 = 118.4*(PI / 180.0);

		bondAngleParams[C_CA_N].k = 63.0;
		bondAngleParams[C_CA_N].l_0 = 110.1*(PI / 180.0);

		bondAngleParams[C_CA_CB].k = 63.0;
		bondAngleParams[C_CA_CB].l_0 = 111.1*(PI / 180.0);

		bondAngleParams[N_CA_CB].k = 80.0;
		bondAngleParams[N_CA_CB].l_0 = 109.5*(PI / 180.0);
		 */
		//ba_ele.batype = POSE::CA_C_N;
		add_bb_angle(ca_0, c_0, n_p1, 116.6*(PI / 180.0), scale_factor * 70.0 );

		//ba_ele.batype = POSE::CA_C_O;
		add_bb_angle(ca_0, c_0, o_0, 120.4*(PI / 180.0), scale_factor * 80.0 );

		//ba_ele.batype = POSE::N_C_O;
		add_bb_angle(n_p1, c_0, o_0,122.9*(PI / 180.0), scale_factor * 80.0 );

		//ba_ele.batype = POSE::C_CA_CB;
		add_bb_angle(c_0, ca_0, cb_0, 111.1*(PI / 180.0), scale_factor * 63.0);

		//ba_ele.batype = POSE::N_CA_CB;
		add_bb_angle(n_0, ca_0, cb_0, 109.5*(PI / 180.0), scale_factor * 80.0);

		//ba_ele.batype = POSE::C_CA_N;
		add_bb_angle(c_0, ca_0, n_0, 110.1*(PI / 180.0), scale_factor * 63.0);

		//ba_ele.batype = POSE::C_N_CA;
		add_bb_angle(c_m1, n_0, ca_0, 121.9*(PI / 180.0), scale_factor * 50.0);

		//ba_ele.batype = POSE::C_N_H;
		add_bb_angle(c_m1, n_0, h_0, 118.4*(PI / 180.0), scale_factor * 38.0);

		//ba_ele.batype = POSE::CA_N_H;
		add_bb_angle(ca_0, n_0, h_0, 119.8*(PI / 180.0) , scale_factor * 35.0);



	}


}




void bb_pose_meta::add_bb_dih(dihedral_element& ele){
	PRODART::POSE::const_atom_shared_ptr at1 = ele.atom1_ptr,
			at2 = ele.atom2_ptr,
			at3 = ele.atom3_ptr,
			at4 = ele.atom4_ptr;
	if (at1 && at2 && at3 && at4){
		ele.atom1 = at1->get_seq_num();
		ele.atom2 = at2->get_seq_num();
		ele.atom3 = at3->get_seq_num();
		ele.atom4 = at4->get_seq_num();
		ele.res_num1 = at1->get_residue()->get_internal_residue_index();
		ele.res_num2 = at2->get_residue()->get_internal_residue_index();
		ele.res_num3 = at3->get_residue()->get_internal_residue_index();
		ele.res_num4 = at4->get_residue()->get_internal_residue_index();
		bb_dih_angle_list.push_back(ele);
	}
	else{

	}
}
void bb_pose_meta::createDihedralList(){

	const double scale_to_joules = 4184; //scale kcal to J
	const double R = 8.314472; //J K-1 mol-1
	const double assumed_temperature = 200.0;
	const double scale_factor = (scale_to_joules / (R * assumed_temperature));
	const double torsion_unit = scale_factor;

	bb_dih_angle_list.resize(0);
	bb_dih_angle_list.reserve(pose_->get_residue_count() * 50);

	POSE::residue_shared_ptr res_0 , res_p1 , res_m1 ;
	POSE::chain_shared_ptr chain_0 , chain_1 , chain_m1 ;
	//cout << "1\n";

	const int resCount = pose_->get_residue_count();//pose_->getResidueCount();

	res_0 = residue_shared_ptr();
	res_m1 = residue_shared_ptr();
	res_m1 = residue_shared_ptr();
	res_p1 = pose_->get_residue(0);
	chain_1 = res_p1->get_chain();
	for (int i = 0; i < resCount ; i++){
		//cout << "2";
		res_m1 = res_0;
		res_0 = res_p1;
		chain_m1 = chain_0;
		chain_0 = chain_1;
		if (i+1 < resCount){
			res_p1 = pose_->get_residue( i + 1 );
		}
		else {
			res_p1 = residue_shared_ptr();
		}

		if (res_p1){
			chain_1 = res_p1->get_chain();
		}
		else {
			chain_1 = chain_shared_ptr();
		}
		POSE::atom_shared_ptr o_m1, ca_m1, ca_0, n_0, h_0, c_m1, o_0, c_0;


		if (res_m1){
			c_m1 = pose_->get_bb_atom(POSE::C, i-1);
			ca_m1 = pose_->get_bb_atom(POSE::CA, i-1);
			o_m1 = pose_->get_bb_atom(POSE::O, i-1);
			//n_m1 = res_m1 pose_->get_bb_atom(POSE::N);
		}
		else {
			c_m1 = atom_shared_ptr();
			ca_m1 = atom_shared_ptr();
			o_m1 = atom_shared_ptr();
			//n_m1 = atom_shared_ptr();
		}


		h_0 = pose_->get_bb_atom(POSE::H, i);
		n_0 = pose_->get_bb_atom(POSE::N, i);
		ca_0 = pose_->get_bb_atom(POSE::CA, i);
		c_0 = pose_->get_bb_atom(POSE::C, i);
		o_0 = pose_->get_bb_atom(POSE::O, i);



		if (chain_0 == chain_m1 && ca_m1  && c_m1  && n_0  && ca_0 ){
			dihedral_element da_ele;
			//da_ele.dihtype = POSE::CA_C_N_CA;
			//da_ele.resNum = i;.
			/*
			tor_ele.amplitude = 2.500 * torsion_unit;
			tor_ele.phase =  180.0 * (PI / 180.0);
			tor_ele.periodicity = 2;
			dihedralAngleParams[CA_C_N_CA].elements.push_back(tor_ele);
			*/
			{
				dih_params_element tor_ele;
				tor_ele.amplitude = 2.500 * torsion_unit;
				tor_ele.phase =  180.0 * (PI / 180.0);
				tor_ele.periodicity = 2;
				da_ele.params.push_back(tor_ele);

			}
			da_ele.atom1_ptr = ca_m1;
			da_ele.atom2_ptr = c_m1;
			da_ele.atom3_ptr = n_0;
			da_ele.atom4_ptr = ca_0;
			add_bb_dih(da_ele);
		}
		if (chain_0 == chain_m1 && ca_m1  && c_m1  && n_0  && h_0 ){
			dihedral_element da_ele;
			//da_ele.dihtype = POSE::CA_C_N_H;
			//da_ele.resNum = i;
			/*
			tor_ele.amplitude = 2.500 * torsion_unit;
			tor_ele.phase =  180.0 * (PI / 180.0);
			tor_ele.periodicity = 2;
			dihedralAngleParams[CA_C_N_H].elements.push_back(tor_ele);
			*/
			{
				dih_params_element tor_ele;
				tor_ele.amplitude = 2.500 * torsion_unit;
				tor_ele.phase =  180.0 * (PI / 180.0);
				tor_ele.periodicity = 2;
				da_ele.params.push_back(tor_ele);

			}
			da_ele.atom1_ptr = ca_m1;
			da_ele.atom2_ptr = c_m1;
			da_ele.atom3_ptr = n_0;
			da_ele.atom4_ptr = h_0;
			add_bb_dih(da_ele);
		}
		if (chain_0 == chain_m1 && o_m1  && c_m1  && n_0  && h_0 ){
			dihedral_element da_ele;
			//da_ele.dihtype = POSE::H_N_C_O; 	//CA_C_N_H; found bug here: wrong dihedral angle, now fixed
			//da_ele.resNum = i;
			/*
				tor_ele.amplitude = 0.650 * torsion_unit;
				tor_ele.phase =  0.0 * (PI / 180.0);
				tor_ele.periodicity = 1;
				dihedralAngleParams[H_N_C_O].elements.push_back(tor_ele);
				tor_ele.amplitude = 2.500 * torsion_unit;
				tor_ele.phase =  180.0 * (PI / 180.0);
				tor_ele.periodicity = 2;
				dihedralAngleParams[H_N_C_O].elements.push_back(tor_ele);
			 */
			{
				dih_params_element tor_ele;
				tor_ele.amplitude = 0.650 * torsion_unit;
				tor_ele.phase =  0.0 * (PI / 180.0);
				tor_ele.periodicity = 1;
				da_ele.params.push_back(tor_ele);
			}
			{
				dih_params_element tor_ele;
				tor_ele.amplitude = 2.500 * torsion_unit;
				tor_ele.phase =  180.0 * (PI / 180.0);
				tor_ele.periodicity = 2;
				da_ele.params.push_back(tor_ele);
			}
			da_ele.atom1_ptr = h_0;
			da_ele.atom2_ptr = n_0;
			da_ele.atom3_ptr = c_m1;
			da_ele.atom4_ptr = o_m1;
			add_bb_dih(da_ele);
		}
		if (chain_0 == chain_m1 && o_m1  && c_m1  && n_0  && ca_0 ){
			dihedral_element da_ele;
			//da_ele.dihtype = POSE::O_C_N_CA;
			//da_ele.resNum = i;
			/*
			tor_ele.amplitude = 2.500 * torsion_unit;
			tor_ele.phase =  180.0 * (PI / 180.0);
			tor_ele.periodicity = 2;
			dihedralAngleParams[O_C_N_CA].elements.push_back(tor_ele);
			*/
			{
				dih_params_element tor_ele;
				tor_ele.amplitude = 2.500 * torsion_unit;
				tor_ele.phase =  180.0 * (PI / 180.0);
				tor_ele.periodicity = 2;
				da_ele.params.push_back(tor_ele);

			}
			da_ele.atom1_ptr = ca_0;
			da_ele.atom2_ptr = n_0;
			da_ele.atom3_ptr = c_m1;
			da_ele.atom4_ptr = o_m1;
			add_bb_dih(da_ele);
		}
		if ( o_0  && c_0  && ca_0  && n_0 ){
			dihedral_element da_ele;
			//da_ele.dihtype = POSE::O_C_CA_N;
			//da_ele.resNum = i;
			/*
			tor_ele.amplitude = 0.100 * torsion_unit; // <- normal value
			tor_ele.phase =   180.0 * (PI / 180.0);
			tor_ele.periodicity = 3;
			dihedralAngleParams[O_C_CA_N].elements.push_back(tor_ele); //affects Psi
			*/
			{
				dih_params_element tor_ele;
				tor_ele.amplitude = 0.100 * torsion_unit; // <- normal value
				tor_ele.phase =   180.0 * (PI / 180.0);
				tor_ele.periodicity = 3;
				da_ele.params.push_back(tor_ele);

			}
			da_ele.atom1_ptr = o_0;
			da_ele.atom2_ptr = c_0;
			da_ele.atom3_ptr = ca_0;
			da_ele.atom4_ptr = n_0;
			add_bb_dih(da_ele);
		}

	}

	//return true;

}


void bb_pose_meta::createBonded14_15Lists(){



	const int atom_count = pose_->get_all_atom_count();

	this->bonded_14_list.resize(0);
	this->bonded_15_list.resize(0);
	this->bonded_14_list.reserve(atom_count * 20);
	this->bonded_15_list.reserve(atom_count * 20);

	for (int i = 0; i < atom_count; i++){
		POSE::atom_shared_ptr at_i = pose_->get_atom(i);
		if (at_i){
			if (at_i->isActive() && at_i->isSet()){
				POSE::residue_shared_ptr res_i = at_i->get_residue();
				POSE::chain_shared_ptr ch_i = at_i->get_chain();
				for (int j = i+1; j < atom_count; j++){
					POSE::atom_shared_ptr at_j = pose_->get_atom(j);
					if (at_j){
						if (at_j->isActive() && at_j->isSet()){
							POSE::residue_shared_ptr res_j = at_j->get_residue();
							POSE::chain_shared_ptr ch_j = at_j->get_chain();

							if (ch_i == ch_j){
								const int seq_sep = std::abs(res_i->get_internal_residue_index() - res_j->get_internal_residue_index());
								const int bond_sep = bb_map->get_bond_sep(seq_sep, at_i->get_type(), at_j->get_type());
								if (bond_sep == 3){
									//1-4
									nb_pair_element ele;
									ele.seq_sep = seq_sep;
									ele.bond_sep = bond_sep;
									ele.atom1_ptr = at_i;
									ele.atom2_ptr = at_j;
									ele.atom1 = ele.atom1_ptr->get_seq_num();
									ele.atom2 = ele.atom2_ptr->get_seq_num();
									ele.atype1 = ele.atom1_ptr->get_type();
									ele.atype2 = ele.atom2_ptr->get_type();
									ele.res_num1 = res_i->get_internal_residue_index();
									ele.res_num2 = res_j->get_internal_residue_index();

									this->bonded_14_list.push_back(ele);

								}
								else if (bond_sep == 4){
									//1-5
									nb_pair_element ele;
									ele.seq_sep = seq_sep;
									ele.bond_sep = bond_sep;
									ele.atom1_ptr = at_i;
									ele.atom2_ptr = at_j;
									ele.atom1 = ele.atom1_ptr->get_seq_num();
									ele.atom2 = ele.atom2_ptr->get_seq_num();
									ele.atype1 = ele.atom1_ptr->get_type();
									ele.atype2 = ele.atom2_ptr->get_type();
									ele.res_num1 = res_i->get_internal_residue_index();
									ele.res_num2 = res_j->get_internal_residue_index();

									this->bonded_15_list.push_back(ele);

								}
							}

						}
					}

				}
			}
		}
	}

}


void bb_pose_meta::update_phi_psi_omega(){

	const int resCount = pose_->get_residue_count();

	for (int i = 0; i < resCount; i++){
		const double phi = pose_->get_phi(i);
		const double psi = pose_->get_psi(i);
		const double omega = pose_->get_omega_to_prev(i);
		double_3_tuple tup(phi, psi, omega);
		this->phi_psi_omega_vec[i] = tup;

		POSE::four_state_sec_struct secs = POSE_UTILS::get_phi_psi_omega_sector(phi, psi, omega);
		residue_conformation_class[i] = secs;
	}

}


void bb_pose_meta::calc_hbond_counts(){
	strand_hbond_count.resize(0);
	strand_hbond_count.resize(pose_->get_residue_count(), 0);
	helix_hbond_count.resize(0);
	helix_hbond_count.resize(pose_->get_residue_count(), 0);

	PRODART::POSE::META::nb_ele_vector& pair_list = this->get_H_O_pair_list();
	PRODART::POSE::four_state_sec_struct_vector& conf_class = this->get_conf_class();
	PRODART::POSE::META::nb_ele_vector::iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		//cout << "\ndb: 1 ";
		if (iter->dist < this->hb_dist_cutoff){
			//cout << "db: 2 ";
			const int seq_sep = iter->seq_sep;
			PRODART::POSE::four_state_sec_struct ss1 = conf_class[iter->res_num1];
			PRODART::POSE::four_state_sec_struct ss2 = conf_class[iter->res_num2];
			//cout << "db: " << iter->atype1.get_label() << " " << iter->atype2.get_label() << " ";
			if ((iter->atype1 == bb_aty_H() &&  iter->atype2 == bb_aty_O())
					|| (iter->atype1 == bb_aty_O() &&  iter->atype2 == bb_aty_H()) ){


				if (seq_sep == 4
						&& (ss1 == PRODART::POSE::ss4_HELIX // secs_HELIX
						&& ss2 == PRODART::POSE::ss4_HELIX)){						// too lax if || not && ?
					helix_hbond_count[iter->res_num1]++;
					helix_hbond_count[iter->res_num2]++;
				}
				else if (seq_sep > 2
						&& (ss1 == PRODART::POSE::ss4_STRAND
						&& ss2 == PRODART::POSE::ss4_STRAND)){						// too lax if || not && ?
					strand_hbond_count[iter->res_num1]++;
					strand_hbond_count[iter->res_num2]++;
				}

			}
		}
	}
}

void bb_pose_meta::calc_sec_struct(){
	const int resCount = pose_->get_residue_count();
	this->residue_sec_struct.resize(0);
	this->residue_sec_struct.resize(resCount, PRODART::POSE::ss3_OTHER);

	for (int i = 1; i < resCount-1; i++){
		if (helix_hbond_count[i] > 0
				&& residue_conformation_class[i] == POSE::ss4_HELIX){
			residue_sec_struct[i] = POSE::ss3_HELIX;
		}
		else if (strand_hbond_count[i] > 0
				&& residue_conformation_class[i] == POSE::ss4_STRAND){
			residue_sec_struct[i] = POSE::ss3_STRAND;
		}
		else if ((residue_conformation_class[i] == POSE::ss4_STRAND
				)
				&& strand_hbond_count[i-1] > 0
				&& strand_hbond_count[i+1] > 0){
			residue_sec_struct[i] = POSE::ss3_STRAND;
		}
		else if ((residue_conformation_class[i] == POSE::ss4_HELIX
				|| helix_hbond_count[i] > 0)
				&& helix_hbond_count[i-1] > 0
				&& helix_hbond_count[i+1] > 0){
			residue_sec_struct[i] = POSE::ss3_HELIX;
		}
		else {
			residue_sec_struct[i] = POSE::ss3_OTHER;//secs_OTHER;
		}
	}

	PRODART::POSE::three_state_sec_struct prev_val = PRODART::POSE::ss3_OTHER;//  secs_OTHER;
	PRODART::POSE::three_state_sec_struct next_val = PRODART::POSE::ss3_OTHER;//  secs_OTHER;
	for (int i = 0; i < resCount - 1; i++){
		next_val = residue_sec_struct[i+1];
		if (residue_sec_struct[i] != POSE::ss3_OTHER){
			if (residue_sec_struct[i] != prev_val
			      && residue_sec_struct[i] != next_val){
				residue_sec_struct[i] = PRODART::POSE::ss3_OTHER;// secs_OTHER;
			}
		}
		prev_val = residue_sec_struct[i];
	}
}

void bb_pose_meta::regularise_atom_names(){

	const int resCount = pose_->get_residue_count();

	/*
	const PRODART::POSE::atom_type bb_aty_N =  PRODART::POSE::atom_type("N");
	const PRODART::POSE::atom_type bb_aty_O =  PRODART::POSE::atom_type("O");
	const PRODART::POSE::atom_type bb_aty_H =  PRODART::POSE::atom_type("H");
	const PRODART::POSE::atom_type bb_aty_CA =  PRODART::POSE::atom_type("CA");
	const PRODART::POSE::atom_type bb_aty_CB =  PRODART::POSE::atom_type("CB");
	const PRODART::POSE::atom_type bb_aty_C =  PRODART::POSE::atom_type("C");
	*/

	for (int i = 0; i < resCount; i++){

		{
			atom_shared_ptr N_at = pose_->get_bb_atom(POSE::N, i);
			if (N_at) N_at->set_type(bb_aty_N());
		}

		{
			atom_shared_ptr O_at = pose_->get_bb_atom(POSE::O, i);
			if (O_at) O_at->set_type(bb_aty_O());
		}
		{
			atom_shared_ptr h_at = pose_->get_bb_atom(POSE::H, i);
			if (h_at) h_at->set_type(bb_aty_H());
		}
		{
			atom_shared_ptr CA_at = pose_->get_bb_atom(POSE::CA, i);
			if (CA_at) CA_at->set_type(bb_aty_CA());
		}
		{
			atom_shared_ptr CB_at = pose_->get_bb_atom(POSE::CB, i);
			if (CB_at) CB_at->set_type(bb_aty_CB());
		}
		{
			atom_shared_ptr C_at = pose_->get_bb_atom(POSE::C, i);
			if (C_at) C_at->set_type(bb_aty_C());
		}
	}

}




}
}
}

