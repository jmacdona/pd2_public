/*
 * bb_strict_omega.cpp
 *
 *  Created on: 4 Nov 2011
 *      Author: jmacdona
 */
#include "bb_strict_omega.h"


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



potential_shared_ptr new_bb_strict_omega(){
	potential_shared_ptr ptr(new bb_strict_omega());
	return ptr;
}



bb_strict_omega::bb_strict_omega(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_strict_omega"));

	b_is_disposable = false;


}

bool bb_strict_omega::init(){

	return true;
}


double bb_strict_omega::get_energy(const double equil_dih_angle, const double dih) const{
	const double energy = std::pow(this->get_dihedral_distance(equil_dih_angle, dih) , 2); //ele.weight *
	return energy;
}

double bb_strict_omega::get_dE_ddih(const double equil_dih_angle, const double dih) const{
	const double dE = 2.0 * (this->get_dihedral_distance(dih, equil_dih_angle)); // ele.weight *
	return dE;
}

double bb_strict_omega::get_energy(POSE::atom_shared_ptr atom1_ptr,
		POSE::atom_shared_ptr atom2_ptr,
		POSE::atom_shared_ptr atom3_ptr,
		POSE::atom_shared_ptr atom4_ptr,
		const double dih,
		const double equil_ang) const{
	if (atom1_ptr
			&& atom2_ptr
			&& atom3_ptr
			&& atom4_ptr){
		if (atom1_ptr->isActive()
				&& atom2_ptr->isActive()
				&& atom3_ptr->isActive()
				&& atom4_ptr->isActive()){
			//const double dih = dihedral(ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords());
			return this->get_energy(equil_ang, dih);
		}
	}

		return 0;

}



double bb_strict_omega::get_energy_ana_grad(POSE::atom_shared_ptr atom1_ptr,
		POSE::atom_shared_ptr atom2_ptr,
		POSE::atom_shared_ptr atom3_ptr,
		POSE::atom_shared_ptr atom4_ptr,
		double dih,
		const double equil_ang,
		UTILS::vector3d_vector& grad, const double weight ) const{
	if (atom1_ptr
			&& atom2_ptr
			&& atom3_ptr
			&& atom4_ptr){
		if (atom1_ptr->isActive()
				&& atom2_ptr->isActive()
				&& atom3_ptr->isActive()
				&& atom4_ptr->isActive()){

			const int indices[4] = {atom1_ptr->get_seq_num(), atom2_ptr->get_seq_num(), atom3_ptr->get_seq_num(), atom4_ptr->get_seq_num()};
			const vector3d init_vec[4] = {atom1_ptr->get_coords(), atom2_ptr->get_coords(), atom3_ptr->get_coords(), atom4_ptr->get_coords()};
			vector3d at_grad_vecs[4], f1;
			//double dih = ele.dih_angle;


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

			const double central_value = get_energy(atom1_ptr, atom2_ptr, atom3_ptr, atom4_ptr, dih, equil_ang);
			const double dE_ddih = weight * this->get_dE_ddih(equil_ang, dih);

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



double bb_strict_omega::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const {
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	for (int i = 0; i < pose_meta_->get_pose()->get_residue_count(); i++){
		POSE::atom_shared_ptr4_tuple tup = pose_meta_->get_pose()->get_omega_atoms(i);
		if (tup.get<0>() && tup.get<1>() && tup.get<2>() && tup.get<3>()){
			if (tup.get<0>()->isActiveAndSet() && tup.get<1>()->isActiveAndSet() && tup.get<2>()->isActiveAndSet() && tup.get<3>()->isActiveAndSet()){
				const double omega = pose_meta_->get_omega(i);
				const double equil_ang = pose_meta_->is_trans(i) ?  -UTILS::PI :  0.0;
				total_energy += this->get_energy(tup.get<0>(), tup.get<1>(), tup.get<2>(), tup.get<3>(), omega, equil_ang);

			}
		}
	}


	/*
	const simple_harmonic_dihedral_element_vector& pair_list = this->dih_vec;//bb_meta_dat->get_bb_dih_angle_list();
	simple_harmonic_dihedral_element_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy(*iter);
	}
	 */

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_strict_omega::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	for (int i = 0; i < pose_meta_->get_pose()->get_residue_count(); i++){
		POSE::atom_shared_ptr4_tuple tup = pose_meta_->get_pose()->get_omega_atoms(i);
		if (tup.get<0>() && tup.get<1>() && tup.get<2>() && tup.get<3>()){
			if (tup.get<0>()->isActiveAndSet() && tup.get<1>()->isActiveAndSet() && tup.get<2>()->isActiveAndSet() && tup.get<3>()->isActiveAndSet()){
				const double omega = pose_meta_->get_omega(i);
				const double equil_ang = pose_meta_->is_trans(i) ?  -UTILS::PI :  0.0;
				total_energy += this->get_energy_ana_grad(tup.get<0>(), tup.get<1>(), tup.get<2>(), tup.get<3>(), omega, equil_ang,
						grad, weight);

			}
		}
	}



	/*
	const simple_harmonic_dihedral_element_vector& pair_list = this->dih_vec;//bb_meta_dat->get_bb_dih_angle_list();
	simple_harmonic_dihedral_element_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		//total_energy += this->get_energy_grad(*iter, grad, weight);
		total_energy += this->get_energy_ana_grad(*iter, grad, weight);
		//this->get_energy_grad(*iter, grad, weight);
		//cout << endl;
	}
	*/

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_strict_omega::get_dihedral_distance(const double dih1, const double dih2)const{

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

/*
void bb_strict_omega::dummy_funct(double val) const{
	//double temp = val;
}
*/





}
}
}
}


