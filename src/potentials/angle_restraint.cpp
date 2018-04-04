/*
 * angle_restraint.cpp
 *
 *  Created on: 21 Feb 2011
 *      Author: jmacdona
 */

#include "angle_restraint.h"

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

potential_shared_ptr new_angle_restraint(){
	potential_shared_ptr ptr(new angle_restraint());
	return ptr;
}

const double angle_restraint::h = 1e-8;

angle_restraint::angle_restraint(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("angle_restraint"));




}

bool angle_restraint::init(){

	return true;
}

double angle_restraint::get_energy(const double ang, const double equilibrium_angle, const double half_angle_const) const{
	return (half_angle_const *
		pow(ang - equilibrium_angle, 2));
}

double angle_restraint::get_energy(const PRODART::POSE::META::angle_element & ele) const{
	if (ele.atom1_ptr->isActiveAndSet()
			&& ele.atom2_ptr->isActiveAndSet()
			&& ele.atom3_ptr->isActiveAndSet()){
		return get_energy(ele.angle, ele.equilibrium_angle, ele.half_angle_const);
	}
	else {
		return 0;
	}
}

double angle_restraint::get_dE_dang(const double ang, const double equilibrium_angle, const double half_angle_const) const{
	return (half_angle_const *
		2.0 * (ang - equilibrium_angle));
}

double angle_restraint::get_energy_ana_grad(const PRODART::POSE::META::angle_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr->isActiveAndSet()
			&& ele.atom2_ptr->isActiveAndSet()
			&& ele.atom3_ptr->isActiveAndSet()){

		const double dE_dang = weight * get_dE_dang(ele.angle, ele.equilibrium_angle, ele.half_angle_const);
		const int indices[3] = {ele.atom1, ele.atom2, ele.atom3};
		const vector3d init_vec[3] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords()};
		//vector3d vec[3] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords()};;


		/*
		cout << ele.atom1_ptr->get_residue()->get_pdb_residue_index() << "\t" << ele.atom1_ptr->get_type().get_label() << "\t"
				<< ele.atom2_ptr->get_residue()->get_pdb_residue_index() << "\t" << ele.atom2_ptr->get_type().get_label() << "\t"
				<< ele.atom3_ptr->get_residue()->get_pdb_residue_index() << "\t" << ele.atom3_ptr->get_type().get_label() << "\t"
				<<  endl;

		cout << ele.atom1_ptr->get_coords() << "\t"
				<< ele.atom2_ptr->get_coords() << "\t"
				<< ele.atom3_ptr->get_coords() << "\t"
				<<  endl;
				*/

		double ang = ele.angle;
		vector3d grad_0, grad_1, grad_2 , f1;

		angle_p1_deriv(init_vec[0],
				init_vec[1],
				init_vec[2],
				ang,
				//f1,
				grad_0);//,
		//f2);
		//cout << ele.angle << "\t" << ang << "\t" << angle(init_vec[0], init_vec[1],init_vec[2]) << endl;
		angle_p2_deriv(init_vec[0],
				init_vec[1],
				init_vec[2],
				ang,
				//f1,
				grad_1);//,
		//f2);
		//cout << ele.angle << "\t" << ang << endl;
		angle_p1_deriv(init_vec[2],
				init_vec[1],
				init_vec[0],
				ang,
				//f1,
				grad_2);//,
		//f2);

		//cout << ele.angle << "\t" << ang << endl;

		grad[indices[0]] += dE_dang * grad_0;
		grad[indices[1]] += dE_dang * grad_1;
		grad[indices[2]] += dE_dang * grad_2;

		//cout << "ana: 0\t" << dE_dang * grad_0 << endl;
		//cout << "ana: 1\t" << dE_dang * grad_1 << endl;
		//cout << "ana: 2\t" << dE_dang * grad_2 << endl;
		const double central_value = get_energy(ang, ele.equilibrium_angle, ele.half_angle_const);

		return central_value;
	}
	else {
		return 0;
	}

}


// TODO replace with analytical gradient calc
double angle_restraint::get_energy_grad(const PRODART::POSE::META::angle_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr->isActiveAndSet()
			&& ele.atom2_ptr->isActiveAndSet()
			&& ele.atom3_ptr->isActiveAndSet()){
		const double central_value = get_energy(ele.angle, ele.equilibrium_angle, ele.half_angle_const);
		const int indices[3] = {ele.atom1, ele.atom2, ele.atom3};
		const vector3d init_vec[3] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords()};
		vector3d vec[3] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords()};;

		for (int v = 0; v <3; v++){
			//cout << "num: " << v << "\t";
			for (int i = 0; i <3; i++){
				double temp = init_vec[v][i] + h;
				dummy_funct(temp);
				const double hh =  temp - vec[v][i];

				vec[v][i] = init_vec[v][i];
				vec[v][i] += hh;
				const double new_angle = angle(vec[0], vec[1], vec[2]);
				const double new_val = this->get_energy(new_angle, ele.equilibrium_angle, ele.half_angle_const);

				vec[v][i] = init_vec[v][i];
				vec[v][i] -= hh;
				const double new_low_angle = angle(vec[0], vec[1], vec[2]);
				const double new_low_val = this->get_energy(new_low_angle, ele.equilibrium_angle, ele.half_angle_const);

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
	else {
		return 0;
	}
}



double angle_restraint::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	//const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;


	const angle_ele_vector& pair_list = pose_meta_->get_angle_harmonic_restraint_list();
	angle_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy(*iter);
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double angle_restraint::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	//const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const angle_ele_vector& pair_list = pose_meta_->get_angle_harmonic_restraint_list();
	angle_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy_ana_grad(*iter, grad, weight);//this->get_energy_grad(*iter, grad, weight);
		//this->get_energy_grad(*iter, grad, weight);
		//cout << endl;
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}



void angle_restraint::dummy_funct(double val) const{
	//double temp = val;
}




}
}
}
