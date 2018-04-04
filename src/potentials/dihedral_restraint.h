/*
 * dihedral_restraint.h
 *
 *  Created on: 25 Jan 2011
 *      Author: jmacdona
 */

#ifndef DIHEDRAL_RESTRAINT_H_
#define DIHEDRAL_RESTRAINT_H_




#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"
#include "utils/dihedral_derivative.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

potential_shared_ptr new_dihedral_restraint();

class dihedral_restraint : public potential_interface {


	friend potential_shared_ptr new_dihedral_restraint();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();


	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

protected:

	static const double h;// = 1e-8;

	dihedral_restraint();
	double get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const;
	double get_dE_ddih(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const;
	double get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele) const;
	double get_energy_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;
	double get_energy_ana_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;





};



inline double dihedral_restraint::get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const{
	const double energy = ele.weight * std::pow(UTILS::get_dihedral_distance(ele.equil_dih_angle, dih) , 2);
	return energy;
}

inline double dihedral_restraint::get_dE_ddih(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, const double dih) const{
	const double dE = ele.weight * 2.0 * (UTILS::get_dihedral_distance(dih, ele.equil_dih_angle));
	return dE;
}

inline double dihedral_restraint::get_energy(const PRODART::POSE::META::simple_harmonic_dihedral_element & ele) const{
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
inline double dihedral_restraint::get_energy_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr
			&& ele.atom2_ptr
			&& ele.atom3_ptr
			&& ele.atom4_ptr){
		if (ele.atom1_ptr->isActive()
				&& ele.atom2_ptr->isActive()
				&& ele.atom3_ptr->isActive()
				&& ele.atom4_ptr->isActive()){

			const int indices[4] = {ele.atom1, ele.atom2, ele.atom3, ele.atom4};
			const UTILS::vector3d init_vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};
			UTILS::vector3d vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};;

			const double central_value = get_energy(ele);

			for (int v = 0; v <4; v++){
				//cout << "num: " << v << "\t";
				for (int i = 0; i <3; i++){
					double temp = init_vec[v][i] + h;
					UTILS::dummy_funct(temp);
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

inline double dihedral_restraint::get_energy_ana_grad(const PRODART::POSE::META::simple_harmonic_dihedral_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.atom1_ptr
			&& ele.atom2_ptr
			&& ele.atom3_ptr
			&& ele.atom4_ptr){
		if (ele.atom1_ptr->isActive()
				&& ele.atom2_ptr->isActive()
				&& ele.atom3_ptr->isActive()
				&& ele.atom4_ptr->isActive()){

			const int indices[4] = {ele.atom1, ele.atom2, ele.atom3, ele.atom4};
			const UTILS::vector3d init_vec[4] = {ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords()};
			UTILS::vector3d at_grad_vecs[4], f1;
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



inline double dihedral_restraint::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;


	const PRODART::POSE::META::simple_harmonic_dihedral_element_vector& pair_list = pose_meta_->get_dihedral_harmonic_restraint_list();//bb_meta_dat->get_bb_dih_angle_list();
	PRODART::POSE::META::simple_harmonic_dihedral_element_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy(*iter);
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

inline double dihedral_restraint::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const PRODART::POSE::META::simple_harmonic_dihedral_element_vector& pair_list = pose_meta_->get_dihedral_harmonic_restraint_list();//bb_meta_dat->get_bb_dih_angle_list();
	PRODART::POSE::META::simple_harmonic_dihedral_element_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		//total_energy += this->get_energy_grad(*iter, grad, weight);
		total_energy += this->get_energy_ana_grad(*iter, grad, weight);
		//this->get_energy_grad(*iter, grad, weight);
		//cout << endl;
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}


}
}
}



#endif /* DIHEDRAL_RESTRAINT_H_ */
