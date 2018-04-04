/*
 * bb_dihedral.cpp
 *
 *  Created on: 27 Oct 2010
 *      Author: jmacdona
 */

#include "bb_dihedral.h"

using namespace boost;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;
using namespace std;

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

potential_shared_ptr new_bb_dihedral(){
	potential_shared_ptr ptr(new bb_dihedral());
	return ptr;
}

const double bb_dihedral::h = 1e-8;

bb_dihedral::bb_dihedral(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_dihedral"));




}

bool bb_dihedral::init(){

	return true;
}


double bb_dihedral::get_energy(const PRODART::POSE::META::dihedral_element& ele, const double dih) const{
	if (ele.atom1_ptr->isActiveAndSet()
			&& ele.atom2_ptr->isActiveAndSet()
			&& ele.atom3_ptr->isActiveAndSet()
			&& ele.atom4_ptr->isActiveAndSet()){
		double energy = 0;
		const dih_params_ele_vector &params = ele.params;
		dih_params_ele_vector::const_iterator it;
		for (it = params.begin(); it != params.end(); it++){
			energy += it->get_energy(dih);
		}
		return energy;
	}
	else {
		return 0;
	}
}
double bb_dihedral::get_energy(const PRODART::POSE::META::dihedral_element & ele) const{
	return this->get_energy(ele, ele.dih_angle);
}

double bb_dihedral::get_dE_ddih(const PRODART::POSE::META::dihedral_element& ele, const double dih) const{
	double energy = 0;
	const dih_params_ele_vector &params = ele.params;
	dih_params_ele_vector::const_iterator it;
	for (it = params.begin(); it != params.end(); it++){
		energy += it->get_dE_ddih(dih);
	}
	return energy;
}

double bb_dihedral::get_dE_ddih(const PRODART::POSE::META::dihedral_element& ele) const{
	return this->get_dE_ddih(ele, ele.dih_angle);
}

// TODO replace with analytical gradient calc
double bb_dihedral::get_energy_grad(const PRODART::POSE::META::dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr->isActiveAndSet()
			&& ele.atom2_ptr->isActiveAndSet()
			&& ele.atom3_ptr->isActiveAndSet()
			&& ele.atom4_ptr->isActiveAndSet()){
		const double central_value = get_energy(ele);
		const int indices[4] = {ele.atom1, ele.atom2, ele.atom3, ele.atom4};
		const vector3d init_vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};
		vector3d vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};;

		for (int v = 0; v <4; v++){
			for (int i = 0; i <3; i++){
				double temp = init_vec[v][i] + h;
				dummy_funct(temp);
				const double hh =  temp - vec[v][i];

				vec[v][i] = init_vec[v][i];
				vec[v][i] += hh;
				const double new_dih = dihedral(vec[0], vec[1], vec[2], vec[3]);
				const double new_val = this->get_energy(ele, new_dih);

				vec[v][i] = init_vec[v][i];
				vec[v][i] -= hh;
				const double new_low_dih = dihedral(vec[0], vec[1], vec[2], vec[3]);
				const double new_low_val = this->get_energy(ele, new_low_dih);

				vec[v][i] = init_vec[v][i];
				const double diff = new_val - new_low_val;
				const double this_grad = diff / (2.0 * hh);
				(grad[indices[v]])[i] += weight * this_grad;
			}
		}

		return central_value;
	}
	else {
		return 0;
	}
}

double bb_dihedral::get_energy_ana_grad(const PRODART::POSE::META::dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr->isActiveAndSet()
			&& ele.atom2_ptr->isActiveAndSet()
			&& ele.atom3_ptr->isActiveAndSet()
			&& ele.atom4_ptr->isActiveAndSet()){

		const int indices[4] = {ele.atom1, ele.atom2, ele.atom3, ele.atom4};
		const vector3d init_vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};

		vector3d at_grad_vecs[4], f1;
		double dih = ele.dih_angle;

		/*
		cout << ele.atom1_ptr->get_residue()->get_pdb_residue_index() << "\t" << ele.atom1_ptr->get_type().get_label() << "\t"
				<< ele.atom2_ptr->get_residue()->get_pdb_residue_index() << "\t" << ele.atom2_ptr->get_type().get_label() << "\t"
				<< ele.atom3_ptr->get_residue()->get_pdb_residue_index() << "\t" << ele.atom3_ptr->get_type().get_label() << "\t"
				<< ele.atom4_ptr->get_residue()->get_pdb_residue_index() << "\t" << ele.atom4_ptr->get_type().get_label() << "\t"
				<<  endl;

		cout << ele.atom1_ptr->get_coords() << "\t"
				<< ele.atom2_ptr->get_coords() << "\t"
				<< ele.atom3_ptr->get_coords() << "\t"
				<< ele.atom4_ptr->get_coords() << "\t"
				<<  endl;
				*/

		//cout << "ENRG: " << central_value << endl << endl;

		dihedral_p1_cosine_deriv(
				init_vec[0],
				init_vec[1],
				init_vec[2],
				init_vec[3],
				dih,
				//f1,
				at_grad_vecs[0]);
		dihedral_p1_cosine_deriv(
				init_vec[3],
				init_vec[2],
				init_vec[1],
				init_vec[0],
				dih,
				//f1,
				at_grad_vecs[3]);
		dihedral_p2_cosine_deriv(
				init_vec[0],
				init_vec[1],
				init_vec[2],
				init_vec[3],
				dih,
				//f1,
				at_grad_vecs[1]);
		dihedral_p2_cosine_deriv(
				init_vec[3],
				init_vec[2],
				init_vec[1],
				init_vec[0],
				dih,
				//f1,
				at_grad_vecs[2]);

		const double dE_ddih = weight * this->get_dE_ddih(ele, dih);

		for (int i = 0; i < 4; i++){
			grad[indices[i]] += dE_ddih * at_grad_vecs[i];
			//cout << "ana: " << i << "\t" << dE_ddih * at_grad_vecs[i] << endl;
		}

		const double central_value = get_energy(ele, dih);

		return central_value;
	}
	else {
		return 0;
	}
}



double bb_dihedral::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;


	const dihedral_ele_vector& pair_list = bb_meta_dat->get_bb_dih_angle_list();
	dihedral_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy(*iter);
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_dihedral::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const dihedral_ele_vector& pair_list = bb_meta_dat->get_bb_dih_angle_list();
	dihedral_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		//total_energy += this->get_energy_grad(*iter, grad, weight);
		total_energy += this->get_energy_ana_grad(*iter, grad, weight);
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}



void bb_dihedral::dummy_funct(double val) const{
	//double temp = val;
}




}
}
}
}

