//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose_meta_interface.h
 *
 *  Created on: 5 Mar 2010
 *      Author: jmacdona
 */

#ifndef POSE_META_INTERFACE_H_
#define POSE_META_INTERFACE_H_

#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include "pose_meta_defs.h"
#include "sim_cell.h"
#include <boost/shared_ptr.hpp>
#include "movers/mover_defs.h"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <map>


namespace PRODART {
namespace POSE {
namespace META {

typedef boost::tuple<double, double> double_2_tuple;
typedef boost::tuple<double, double, double> double_3_tuple;
typedef std::vector<double_3_tuple> double_3_tuple_vector;

class pose_meta_interface;

typedef boost::shared_ptr<pose_meta_interface> pose_meta_shared_ptr;
typedef boost::shared_ptr<const pose_meta_interface> const_pose_meta_shared_ptr;
typedef boost::tuple<int, int> int_int_tuple;
typedef boost::tuple<bool, double> bool_double_tuple;
typedef std::map<int_int_tuple, bool_double_tuple> int_int_tup_bool_double_tup_map;

//! class for storing pose meta data (pair lists etc) for use during simulations
class pose_meta_interface {



protected:
	pose_meta_interface(PRODART::POSE::pose_shared_ptr pose_to_attach);
	virtual ~pose_meta_interface(){
	}

	PRODART::POSE::pose_shared_ptr pose_;
	bool refresh_pair_lists;
	double move_margin;

	std::map<POSE::atom_shared_ptr, PRODART::UTILS::vector3d> update_state_store;
	std::map<POSE::atom_shared_ptr, PRODART::UTILS::vector3d> has_moved_coord_store;
	void make_update_state_store();
	void reset_update_state_store();
	bool pair_lists_ok();

	void make_has_moved_coord_store();
	void reset_has_moved_coord_store();
	void update_has_moved_flags();

	bonded_ele_vector bond_list;
	simple_harmonic_dihedral_element_vector dih_restraint_list;
	angle_ele_vector ang_restraint_list;
	coord_rst_element_vector coord_restraint_list;
	sse_axis_element_vector sse_axis_restraint_list;
	upper_lower_bonded_pair_element_vector residue_contact_restraint_list;
	upper_lower_bonded_pair_element_vector ca_GO_restraint_list;
	bool go_restraints_set;


	PRODART::POSE::four_state_sec_struct_vector sec_struct_restraint_list;
	double_vector sec_struct_restraint_weight_list;

	int_int_tup_bool_double_tup_map strand_residue_pairing_rst_map;

	PRODART::POSE::MOVERS::mover_flags mover_return_flags;

	PRODART::UTILS::vector3d_vector gradient;

	PRODART::POSE::four_state_sec_struct_vector residue_conformation_class;
	PRODART::POSE::three_state_sec_struct_vector residue_sec_struct;
	double_3_tuple_vector phi_psi_omega_vec;

	int_vector strand_hbond_count, helix_hbond_count;

	//! update pair list
	virtual void update_pair_lists() = 0;

public:

	enum poseMetaType{pmInterface,
					pm_bb_pose_meta,
					pm_ca_pose_meta};

	virtual poseMetaType get_pose_meta_type() const{
		return pmInterface;
	}

	//! recalc pairs in pair list
	virtual void recalc_pair_lists_dists() = 0;

	//! flag pair list for update on next refresh
	void set_update_pair_lists_flag(){
		this->refresh_pair_lists = true;
	}

	PRODART::POSE::pose_shared_ptr get_pose();

	bonded_ele_vector& get_dist_harmonic_restraint_list(){ //get_bond_list
		return bond_list;
	}

	PRODART::POSE::MOVERS::mover_flags& get_mover_flags(){
		return this->mover_return_flags;
	}

	const PRODART::UTILS::vector3d_vector& get_gradient() const{
		return gradient;
	}

	PRODART::UTILS::vector3d_vector& get_gradient(){
		return gradient;
	}


	PRODART::POSE::four_state_sec_struct_vector& get_conf_class(){
		return residue_conformation_class;
	}

	PRODART::POSE::four_state_sec_struct get_residue_conf_class(const int res_num) const{
		return residue_conformation_class[res_num];
	}

	void set_residue_conf_class(const int res_num, const PRODART::POSE::four_state_sec_struct secstruct_class){
		residue_conformation_class[res_num] = secstruct_class;
	}

	PRODART::POSE::three_state_sec_struct_vector& get_sec_struct(){
		return this->residue_sec_struct;
	}

	PRODART::POSE::three_state_sec_struct get_residue_sec_struct(const int res_num) const{
		return residue_sec_struct[res_num];
	}

	void set_residue_sec_struct(const int res_num, const PRODART::POSE::three_state_sec_struct secstruct){
		residue_sec_struct[res_num] = secstruct;
	}

	double_3_tuple get_phi_psi_omega(const int res_num) const{
		return this->phi_psi_omega_vec[res_num];
	}

	void set_phi_psi_omega(const int res_num, const double_3_tuple tup ){
		phi_psi_omega_vec[res_num] = tup;
	}

	double get_phi(const int res_num) const{
		return this->phi_psi_omega_vec[res_num].get<0>();
	}

	double get_psi(const int res_num) const{
		return this->phi_psi_omega_vec[res_num].get<1>();
	}

	static int get_dihedral_grid_ref(const double dih, const int griddivisions ){
		const double increment = (2.0*UTILS::PI) / (double)griddivisions;
		return static_cast<int>((UTILS::PI + dih) / increment);
	}

	int get_phi_grid_ref( const int res_num, const int griddivisions ) const{

		const double increment = (2*UTILS::PI) / griddivisions;
		const double phi = UTILS::PI + get_phi(res_num);
		return static_cast<int>(phi / increment);

	}

	int get_psi_grid_ref( const int res_num, const int griddivisions ) const{

		const double increment = (2*UTILS::PI) / griddivisions;
		const double psi = UTILS::PI + get_psi(res_num);
		return static_cast<int>(psi / increment);

	}

	static int get_phi_psi_sector(const double phi, const double psi, const int griddivisions ){
		return get_dihedral_grid_ref( phi, griddivisions ) + (griddivisions * get_dihedral_grid_ref( psi, griddivisions ));
	}

	int get_phi_psi_sector( const int res_num, const int griddivisions ) const{
		return get_phi_grid_ref( res_num, griddivisions ) + (griddivisions * get_psi_grid_ref( res_num, griddivisions ));
	}


	double get_omega(const int res_num) const{
		return this->phi_psi_omega_vec[res_num].get<2>();
	}

	bool is_trans(const int res_num) const {
		const double omega = this->get_omega(res_num);
		if ( omega < (UTILS::PI / 2.0) && omega > -(UTILS::PI / 2.0)) {
			//is cis
			return false;
		}
		else{
			//is trans
			return true;
		}
	}

	void set_phi(const int res_num, const double val) {
		this->phi_psi_omega_vec[res_num].get<0>() = val;
	}

	void set_psi(const int res_num, const double val) {
		this->phi_psi_omega_vec[res_num].get<1>() = val;
	}

	void set_omega(const int res_num, const double val) {
		this->phi_psi_omega_vec[res_num].get<2>() = val;
	}

	void clear_dist_harmonic_restraints(){
		this->bond_list.resize(0);
	}

	void add_dist_harmonic_restraint(const int at1_seq_num,
			const int at2_seq_num,
			const double eq_dist,
			const double half_bond_const){

		const int atom_count = pose_->get_all_atom_count();

		if (at1_seq_num < atom_count && at2_seq_num < atom_count){
			POSE::atom_shared_ptr at1 = this->pose_->get_atom(at1_seq_num);
			POSE::atom_shared_ptr at2 = this->pose_->get_atom(at2_seq_num);

			if (at1 && at2){
				bonded_pair_element ele;
				ele.atom1_ptr = at1;
				ele.atom2_ptr = at2;
				ele.atom1 = at1->get_seq_num();
				ele.atom2 = at2->get_seq_num();
				ele.equilibrium_dist = eq_dist;
				ele.half_bond_const = half_bond_const;
				ele.dist_sq = (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()).mod_sq();
				ele.dist = std::sqrt(ele.dist_sq);
				this->bond_list.push_back(ele);


				if ((!at1->isActive()) || (!at2->isActive())){
					std::cerr << "WARNING: pose_meta_interface: harmonic_restraint: atoms not active so will have no effect" << std::endl;
				}

			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add harmonic_restraint: pointers are NULL" << std::endl;
			}
		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add harmonic_restraint: sequence numbers not valid" << std::endl;
		}

	}

	void clear_ca_GO_restraints(){
		const int residue_count = pose_->get_residue_count();
		this->ca_GO_restraint_list.resize(0);
		upper_lower_bonded_pair_element ele;
		ele.is_valid = false;
		this->ca_GO_restraint_list.resize((residue_count*residue_count), ele);
		go_restraints_set = false;
	}

	/*
	upper_lower_bonded_pair_element_vector& get_ca_GO_restraint_list(){
		return this->ca_GO_restraint_list;
	}
	*/
	upper_lower_bonded_pair_element get_ca_GO_restraint(const int res1_num,
			const int res2_num){
		const int residue_count = pose_->get_residue_count();

		if (res1_num < residue_count && res2_num < residue_count){
			const int res1 = std::min(res1_num, res2_num);
			const int res2 = std::max(res1_num, res2_num);
			const unsigned int index = res1 + (res2*(pose_->get_residue_count()));
			if (index < this->ca_GO_restraint_list.size()){
				return ca_GO_restraint_list[index];
			}
			else {
				std::cerr << "WARNING: pose_meta_interface: add_ca_GO_restraint: index is out of range" << std::endl;
			}

		}
		else {
			std::cerr << "ERROR: pose_meta_interface: get_ca_GO_restraint_list: sequence numbers not valid" << std::endl;

		}

		upper_lower_bonded_pair_element ele;
		ele.is_valid = false;
		return ele;


	}

	void add_ca_GO_restraint(int res1_num,
			int res2_num,
			upper_lower_bonded_pair_element& ele){
		const int residue_count = pose_->get_residue_count();

		if (res1_num > res2_num){
			std::swap(res1_num,res2_num);
		}

		if (res1_num < residue_count && res2_num < residue_count){
			POSE::atom_shared_ptr at1 = this->pose_->get_bb_atom(POSE::CA,res1_num);
			POSE::atom_shared_ptr at2 = this->pose_->get_bb_atom(POSE::CA, res2_num);

			if (at1 && at2){
				ele.res_num1 = res1_num;
				ele.res_num2 = res2_num;
				ele.atom1_ptr = at1;
				ele.atom2_ptr = at2;
				ele.atom1 = at1->get_seq_num();
				ele.atom2 = at2->get_seq_num();
				ele.dist_sq = (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()).mod_sq();
				ele.dist = std::sqrt(ele.dist_sq);
				ele.is_valid = true;
				const unsigned int index = res1_num + (res2_num*(pose_->get_residue_count()));
				if (index < this->ca_GO_restraint_list.size()){
					this->ca_GO_restraint_list[index] = ele;
					go_restraints_set = true;

					//std::cerr << "INFO: pose_meta_interface: add_ca_GO_restraint: add rst" << std::endl;

				}
				else {
					std::cerr << "WARNING: pose_meta_interface: add_ca_GO_restraint: index is out of range" << std::endl;
				}

				if ((!at1->isActive()) || (!at2->isActive())){
					std::cerr << "WARNING: pose_meta_interface: add_ca_GO_restraint: atoms not active so will have no effect" << std::endl;
				}

			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add add_ca_GO_restraint: pointers are NULL" << std::endl;
			}

		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add add_ca_GO_restraint: sequence numbers not valid" << std::endl;
		}
	}

	void add_ca_GO_restraint(int res1_num,
			int res2_num,
			const double eq_dist_lower,
			const double eq_dist_upper,
			const double half_bond_const){

		//const int atom_count = pose_->get_all_atom_count();
		const int residue_count = pose_->get_residue_count();

		if (res1_num > res2_num){
			std::swap(res1_num,res2_num);
		}

		if (res1_num < residue_count && res2_num < residue_count){
			POSE::atom_shared_ptr at1 = this->pose_->get_bb_atom(POSE::CA,res1_num);
			POSE::atom_shared_ptr at2 = this->pose_->get_bb_atom(POSE::CA, res2_num);

			if (at1 && at2){
				upper_lower_bonded_pair_element ele;
				ele.res_num1 = res1_num;
				ele.res_num2 = res2_num;
				ele.atom1_ptr = at1;
				ele.atom2_ptr = at2;
				ele.atom1 = at1->get_seq_num();
				ele.atom2 = at2->get_seq_num();
				ele.equilibrium_dist_lower = eq_dist_lower;
				ele.equilibrium_dist_upper = eq_dist_upper;
				ele.half_bond_const = half_bond_const;
				ele.dist_sq = (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()).mod_sq();
				ele.dist = std::sqrt(ele.dist_sq);
				ele.is_valid = true;
				const unsigned int index = res1_num + (res2_num*(pose_->get_residue_count()));
				if (index < this->ca_GO_restraint_list.size()){
					this->ca_GO_restraint_list[index] = ele;
					go_restraints_set = true;

					//std::cerr << "INFO: pose_meta_interface: add_ca_GO_restraint: add rst" << std::endl;

				}
				else {
					std::cerr << "WARNING: pose_meta_interface: add_ca_GO_restraint: index is out of range" << std::endl;
				}

				if ((!at1->isActive()) || (!at2->isActive())){
					std::cerr << "WARNING: pose_meta_interface: add_ca_GO_restraint: atoms not active so will have no effect" << std::endl;
				}

			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add add_ca_GO_restraint: pointers are NULL" << std::endl;
			}
		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add add_ca_GO_restraint: sequence numbers not valid" << std::endl;
		}

	}

	bool are_go_restraints_set() const{
		return go_restraints_set;
	}


	void clear_residue_contact_restraints(){
		this->residue_contact_restraint_list.resize(0);
	}

	upper_lower_bonded_pair_element_vector& get_residue_contact_restraint_list(){
		return this->residue_contact_restraint_list;
	}


	void add_residue_contact_restraint(const int res1_num,
			const int res2_num,
			const double eq_dist_lower,
			const double eq_dist_upper,
			const double half_bond_const){

		//const int atom_count = pose_->get_all_atom_count();
		const int residue_count = pose_->get_residue_count();

		if (res1_num < residue_count && res2_num < residue_count){
			POSE::atom_shared_ptr at1 = this->pose_->get_bb_atom(POSE::CA,res1_num);
			POSE::atom_shared_ptr at2 = this->pose_->get_bb_atom(POSE::CA, res2_num);

			if (at1 && at2){
				upper_lower_bonded_pair_element ele;
				ele.res_num1 = res1_num;
				ele.res_num2 = res2_num;
				ele.atom1_ptr = at1;
				ele.atom2_ptr = at2;
				ele.atom1 = at1->get_seq_num();
				ele.atom2 = at2->get_seq_num();
				ele.equilibrium_dist_lower = eq_dist_lower;
				ele.equilibrium_dist_upper = eq_dist_upper;
				ele.half_bond_const = half_bond_const;
				ele.dist_sq = (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()).mod_sq();
				ele.dist = std::sqrt(ele.dist_sq);
				this->residue_contact_restraint_list.push_back(ele);


				if ((!at1->isActive()) || (!at2->isActive())){
					std::cerr << "WARNING: pose_meta_interface: harmonic_restraint: atoms not active so will have no effect" << std::endl;
				}

			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add harmonic_restraint: pointers are NULL" << std::endl;
			}
		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add harmonic_restraint: sequence numbers not valid" << std::endl;
		}

	}

	void add_sse_restraint(const int start_resnum,
			const int end_resnum,
			POSE::four_state_sec_struct sse_type,
			const double half_bond_const,
			const UTILS::vector3d start,
			const UTILS::vector3d end);

	void clear_sse_axis_restraint_list(){
		this->sse_axis_restraint_list.clear();
	}

	sse_axis_element_vector& get_sse_axis_restraint_list(){
		return this->sse_axis_restraint_list;
	}

	PRODART::POSE::four_state_sec_struct_vector& get_sec_struct_restraint_list(){
		return sec_struct_restraint_list;
	}

	double_vector& get_sec_struct_restraint_weight_list(){
		return sec_struct_restraint_weight_list;
	}

	void clear_sec_struct_restraint_list(){
		sec_struct_restraint_list.resize(0);
		sec_struct_restraint_list.resize(pose_->get_residue_count(), POSE::ss4_UNDEF);
		sec_struct_restraint_weight_list.resize(0);
		sec_struct_restraint_weight_list.resize(pose_->get_residue_count(),0.0);
	}

	void add_sec_struct_restraint(const int resnum, const POSE::four_state_sec_struct val, const double weight){
		if (resnum >= pose_->get_residue_count()){
			std::cerr << "ERROR: pose_meta_interface: can't add sec_struct_restraint: residue number is out of range" << std::endl;
		}
		else {
			sec_struct_restraint_list[resnum] = val;
			sec_struct_restraint_weight_list[resnum] = weight;
		}
	}

	simple_harmonic_dihedral_element_vector& get_dihedral_harmonic_restraint_list(){
		return dih_restraint_list;
	}

	void clear_dihedral_harmonic_restraints(){
		this->dih_restraint_list.resize(0);
	}

	void add_strand_residue_pair_restraint(const int resnum1, const int resnum2, const bool is_parallel, const double weight){

		const int rescount = pose_->get_residue_count();

		if (resnum1 < rescount && resnum2 < rescount){

			int_int_tuple res_tup(resnum1, resnum2);
			bool_double_tuple spec_tup(is_parallel, weight);

			strand_residue_pairing_rst_map[res_tup] = spec_tup;

		}
		else {
			std::cerr << "ERROR: pose_meta_interface: can't add strand_residue_pair_restraint: residue numbers not valid" << std::endl;
		}



	}

	int_int_tup_bool_double_tup_map& get_strand_residue_pair_restraints(){
		return strand_residue_pairing_rst_map;
	}

	void clear_strand_residue_pair_restraints(){
		this->strand_residue_pairing_rst_map.clear();
	}

	void add_dihedral_harmonic_restraint(const int at1_seq_num,
			const int at2_seq_num,
			const int at3_seq_num,
			const int at4_seq_num,
			const double eq_dih,
			const double weight){
		const int atom_count = pose_->get_all_atom_count();

		if (at1_seq_num < atom_count
				&& at2_seq_num < atom_count
				&& at3_seq_num <atom_count
				&& at4_seq_num < atom_count){
			POSE::atom_shared_ptr at1 = this->pose_->get_atom(at1_seq_num);
			POSE::atom_shared_ptr at2 = this->pose_->get_atom(at2_seq_num);
			POSE::atom_shared_ptr at3 = this->pose_->get_atom(at3_seq_num);
			POSE::atom_shared_ptr at4 = this->pose_->get_atom(at4_seq_num);

			if (at1 && at2 && at3 && at4){
				simple_harmonic_dihedral_element ele;
				ele.atom1_ptr = at1;
				ele.atom2_ptr = at2;
				ele.atom3_ptr = at3;
				ele.atom4_ptr = at4;
				ele.atom1 = at1->get_seq_num();
				ele.atom2 = at2->get_seq_num();
				ele.atom3 = at3->get_seq_num();
				ele.atom4 = at4->get_seq_num();
				ele.equil_dih_angle = eq_dih;
				ele.weight = weight;
				ele.dih_angle = UTILS::dihedral(ele.atom1_ptr->get_coords(),
						ele.atom2_ptr->get_coords(),
						ele.atom3_ptr->get_coords(),
						ele.atom4_ptr->get_coords());
				this->dih_restraint_list.push_back(ele);


				if ((!at1->isActive())
						|| (!at2->isActive())
						|| (!at3->isActive())
						|| (!at4->isActive())){
					std::cerr << "WARNING: pose_meta_interface: add_dihedral_harmonic_restraint: atoms not active so will have no effect" << std::endl;
				}

			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add dihedral_harmonic_restraint: pointers are NULL" << std::endl;
			}
		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add dihedral_harmonic_restraint: sequence numbers not valid" << std::endl;
		}
	}

	void clear_angle_harmonic_restraints(){
		ang_restraint_list.resize(0);
	}

	void add_angle_harmonic_restraint(const int at1_seq_num,
			const int at2_seq_num,
			const int at3_seq_num,
			const double eq_dih,
			const double weight){
		const int atom_count = pose_->get_all_atom_count();

		if (at1_seq_num < atom_count
				&& at2_seq_num < atom_count
				&& at3_seq_num <atom_count){
			POSE::atom_shared_ptr at1 = this->pose_->get_atom(at1_seq_num);
			POSE::atom_shared_ptr at2 = this->pose_->get_atom(at2_seq_num);
			POSE::atom_shared_ptr at3 = this->pose_->get_atom(at3_seq_num);


			if (at1 && at2 && at3 ){
				angle_element ele;
				ele.atom1_ptr = at1;
				ele.atom2_ptr = at2;
				ele.atom3_ptr = at3;

				ele.atom1 = at1->get_seq_num();
				ele.atom2 = at2->get_seq_num();
				ele.atom3 = at3->get_seq_num();

				ele.equilibrium_angle = eq_dih;
				ele.half_angle_const = weight;
				ele.angle = UTILS::angle(ele.atom1_ptr->get_coords(),
						ele.atom2_ptr->get_coords(),
						ele.atom3_ptr->get_coords());
				this->ang_restraint_list.push_back(ele);


				if ((!at1->isActive())
						|| (!at2->isActive())
						|| (!at3->isActive())){
					std::cerr << "WARNING: pose_meta_interface: add_dihedral_harmonic_restraint: atoms not active so will have no effect" << std::endl;
				}

			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add dihedral_harmonic_restraint: pointers are NULL" << std::endl;
			}
		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add dihedral_harmonic_restraint: sequence numbers not valid" << std::endl;
		}
	}

	angle_ele_vector& get_angle_harmonic_restraint_list(){
		return ang_restraint_list;
	}

	void add_harmonic_coord_restraint(const int at_seq_num,
			const UTILS::vector3d pos,
			const double half_bond_const,
			const double eq_dist = 0){
		const int atom_count = pose_->get_all_atom_count();

		if (at_seq_num < atom_count ){
			POSE::atom_shared_ptr at = this->pose_->get_atom(at_seq_num);

			if (at){
				coord_rst_element ele;
				ele.atom_ptr = at;
				ele.atom = at->get_seq_num();
				ele.coord = pos;
				ele.equilibrium_dist = eq_dist;
				ele.half_bond_const = half_bond_const;
				ele.dist_sq = (ele.atom_ptr->get_coords() - ele.coord).mod_sq();
				ele.dist = std::sqrt(ele.dist_sq);
				this->coord_restraint_list.push_back(ele);


				if ((!at->isActive())){
					std::cerr << "WARNING: pose_meta_interface: add_harmonic_coord_restraint: atom not active so will have no effect" << std::endl;
				}

			}
			else{
				std::cerr << "ERROR: pose_meta_interface: can't add_harmonic_coord_restraint: pointer is NULL" << std::endl;
			}
		}
		else{
			std::cerr << "ERROR: pose_meta_interface: can't add_harmonic_coord_restraint: sequence numbers not valid" << std::endl;
		}
	}

	coord_rst_element_vector& get_coord_restraint_list() {
		return coord_restraint_list;
	}


protected:

	void updateBondedPairList(bonded_ele_vector& bond_list);
	void updateAngleList(angle_ele_vector& angle_list);
	void updateDihedralList(dihedral_ele_vector& dih_list);
	void updateDihedralList(simple_harmonic_dihedral_element_vector& dih_list);
	void updateNonBondedPairList(nb_ele_vector& bond_list);
	void updateCoordRestraintList(coord_rst_element_vector& coord_rst_list);

private:

};

inline PRODART::POSE::pose_shared_ptr pose_meta_interface::get_pose(){
	return pose_;
}




}
}
}



#endif /* POSE_META_INTERFACE_H_ */
