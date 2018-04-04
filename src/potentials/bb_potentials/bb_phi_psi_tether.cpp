/*
 * bb_phi_psi_tether.cpp
 *
 *  Created on: 29 Oct 2010
 *      Author: jmacdona
 */
#include "bb_phi_psi_tether.h"

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

potential_shared_ptr new_bb_phi_psi_tether(){
	potential_shared_ptr ptr(new bb_phi_psi_tether());
	return ptr;
}

const double bb_phi_psi_tether::h = 1e-8;

bb_phi_psi_tether::bb_phi_psi_tether(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_phi_psi_tether"));
	equil_set = false;
	b_is_disposable = true;


}

bool bb_phi_psi_tether::init(){

	return true;
}


double bb_phi_psi_tether::get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const{
	const double energy = ele.weight * std::pow(this->get_dihedral_distance(ele.equil_dih_angle, dih) , 2);
	return energy;
}

double bb_phi_psi_tether::get_dE_ddih(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const{
	const double dE = ele.weight * 2.0 * (this->get_dihedral_distance(dih, ele.equil_dih_angle));
	return dE;
}

double bb_phi_psi_tether::get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element & ele) const{
	if (ele.atom1_ptr
			&& ele.atom2_ptr
			&& ele.atom3_ptr
			&& ele.atom4_ptr){
		if (ele.atom1_ptr->isActive()
				&& ele.atom2_ptr->isActive()
				&& ele.atom3_ptr->isActive()
				&& ele.atom4_ptr->isActive()){
			const double dih = dihedral(ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords());
			return this->get_energy(ele, dih);
		}
	}

		return 0;

}

// TODO replace with analytical gradient calc
double bb_phi_psi_tether::get_energy_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr
			&& ele.atom2_ptr
			&& ele.atom3_ptr
			&& ele.atom4_ptr){
		if (ele.atom1_ptr->isActive()
				&& ele.atom2_ptr->isActive()
				&& ele.atom3_ptr->isActive()
				&& ele.atom4_ptr->isActive()){

			const int indices[4] = {ele.atom1, ele.atom2, ele.atom3, ele.atom4};
			const vector3d init_vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};
			vector3d vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};;

			const double central_value = get_energy(ele);

			for (int v = 0; v <4; v++){
				//cout << "num: " << v << "\t";
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
					//cout << weight * this_grad << "\t";
				}
				//cout << endl;
			}

			return central_value;
		}
	}

	return 0;

}

double bb_phi_psi_tether::get_energy_ana_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr
			&& ele.atom2_ptr
			&& ele.atom3_ptr
			&& ele.atom4_ptr){
		if (ele.atom1_ptr->isActive()
				&& ele.atom2_ptr->isActive()
				&& ele.atom3_ptr->isActive()
				&& ele.atom4_ptr->isActive()){

			const int indices[4] = {ele.atom1, ele.atom2, ele.atom3, ele.atom4};
			const vector3d init_vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};
			vector3d at_grad_vecs[4], f1;
			double dih = ele.dih_angle;


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

			const double central_value = get_energy(ele, dih);
			const double dE_ddih = weight * this->get_dE_ddih(ele, dih);

			/*
			cout << radians_to_degrees(ele.dih_angle) << "\t"
					<< radians_to_degrees(dihedral(init_vec[0],init_vec[1],init_vec[2],init_vec[3])) << "\t"
					<< radians_to_degrees(dih) << "\t:\t"
					<< radians_to_degrees(ele.equil_dih_angle) << "\t"
					<< radians_to_degrees(this->get_dihedral_distance(dih, ele.equil_dih_angle)) << endl;
					*/

			for (int i = 0; i < 4; i++){
				grad[indices[i]] += dE_ddih * at_grad_vecs[i];
				//cout << "ana: " << i << "\t" << dE_ddih * at_grad_vecs[i] << endl;
			}


			return central_value;
		}
	}

	return 0;
}



double bb_phi_psi_tether::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	if (!this->equil_set) this->set_equil(pose_meta_);
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;


	const simple_harmonic_dihedral_element_vector& pair_list = this->dih_vec;//bb_meta_dat->get_bb_dih_angle_list();
	simple_harmonic_dihedral_element_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy(*iter);
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_phi_psi_tether::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	if (!this->equil_set) this->set_equil(pose_meta_);
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const simple_harmonic_dihedral_element_vector& pair_list = this->dih_vec;//bb_meta_dat->get_bb_dih_angle_list();
	simple_harmonic_dihedral_element_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		//total_energy += this->get_energy_grad(*iter, grad, weight);
		total_energy += this->get_energy_ana_grad(*iter, grad, weight);
		//this->get_energy_grad(*iter, grad, weight);
		//cout << endl;
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_phi_psi_tether::get_dihedral_distance(const double dih1, const double dih2)const{

	const double diff = dih1 - dih2;

	if (fabs(diff) < PI){
		return diff;
	}
	else if (diff > 0){
		return (diff - (2.0 * PI));
	}
	else {
		return (diff + (2.0 * PI));
	}

}


void bb_phi_psi_tether::dummy_funct(double val) const{
	//double temp = val;
}

void bb_phi_psi_tether::set_equil(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_) const{
	equil_set = true;
	equil_vals.clear();
	dih_vec.clear();
	equil_vals.reserve(pose_meta_->get_pose()->get_residue_count());
	dih_vec.reserve(pose_meta_->get_pose()->get_residue_count());
	for (int i = 0; i < pose_meta_->get_pose()->get_residue_count(); i++){
		equil_vals.push_back(pose_meta_->get_phi_psi_omega(i));
		{
			META::simple_harmonic_dihedral_element ele;
			POSE::atom_shared_ptr4_tuple tup = pose_meta_->get_pose()->get_phi_atoms(i);
			ele.atom1_ptr = tup.get<0>();
			ele.atom2_ptr = tup.get<1>();
			ele.atom3_ptr = tup.get<2>();
			ele.atom4_ptr = tup.get<3>();
			if (ele.atom1_ptr
					&& ele.atom2_ptr
					&& ele.atom3_ptr
					&& ele.atom4_ptr){
				if (ele.atom1_ptr->isActive()
						&& ele.atom2_ptr->isActive()
						&& ele.atom3_ptr->isActive()
						&& ele.atom4_ptr->isActive()){
					ele.atom1 = ele.atom1_ptr->get_seq_num();
					ele.atom2 = ele.atom2_ptr->get_seq_num();
					ele.atom3 = ele.atom3_ptr->get_seq_num();
					ele.atom4 = ele.atom4_ptr->get_seq_num();
				}
			}
			ele.equil_dih_angle = equil_vals[i].get<0>();
			ele.dih_angle = ele.equil_dih_angle;
			ele.weight = 1.0;
			dih_vec.push_back(ele);
		}
		{
			META::simple_harmonic_dihedral_element ele;
			POSE::atom_shared_ptr4_tuple tup = pose_meta_->get_pose()->get_psi_atoms(i);
			ele.atom1_ptr = tup.get<0>();
			ele.atom2_ptr = tup.get<1>();
			ele.atom3_ptr = tup.get<2>();
			ele.atom4_ptr = tup.get<3>();
			if (ele.atom1_ptr
					&& ele.atom2_ptr
					&& ele.atom3_ptr
					&& ele.atom4_ptr){
				if (ele.atom1_ptr->isActive()
						&& ele.atom2_ptr->isActive()
						&& ele.atom3_ptr->isActive()
						&& ele.atom4_ptr->isActive()){
					ele.atom1 = ele.atom1_ptr->get_seq_num();
					ele.atom2 = ele.atom2_ptr->get_seq_num();
					ele.atom3 = ele.atom3_ptr->get_seq_num();
					ele.atom4 = ele.atom4_ptr->get_seq_num();
				}
			}
			ele.equil_dih_angle = equil_vals[i].get<1>();
			ele.dih_angle = ele.equil_dih_angle;
			ele.weight = 1.0;
			dih_vec.push_back(ele);
		}
		{
			META::simple_harmonic_dihedral_element ele;
			POSE::atom_shared_ptr4_tuple tup = pose_meta_->get_pose()->get_omega_atoms(i);
			ele.atom1_ptr = tup.get<0>();
			ele.atom2_ptr = tup.get<1>();
			ele.atom3_ptr = tup.get<2>();
			ele.atom4_ptr = tup.get<3>();
			if (ele.atom1_ptr
					&& ele.atom2_ptr
					&& ele.atom3_ptr
					&& ele.atom4_ptr){
				if (ele.atom1_ptr->isActive()
						&& ele.atom2_ptr->isActive()
						&& ele.atom3_ptr->isActive()
						&& ele.atom4_ptr->isActive()){
					ele.atom1 = ele.atom1_ptr->get_seq_num();
					ele.atom2 = ele.atom2_ptr->get_seq_num();
					ele.atom3 = ele.atom3_ptr->get_seq_num();
					ele.atom4 = ele.atom4_ptr->get_seq_num();
				}
			}
			ele.equil_dih_angle = equil_vals[i].get<2>();
			ele.dih_angle = ele.equil_dih_angle;
			if (pose_meta_->is_trans(i)){
				ele.equil_dih_angle = -UTILS::PI;
			}
			else{
				ele.equil_dih_angle = 0;
			}

			ele.weight = 30.0;
			dih_vec.push_back(ele);
		}

	}

}


}
}
}
}

