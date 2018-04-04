/*
 * lj_nb_potential.h
 *
 *  Created on: 22 Nov 2010
 *      Author: jmacdona
 */

#ifndef LJ_NB_POTENTIAL_H_
#define LJ_NB_POTENTIAL_H_


#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <map>

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

typedef boost::tuple<double, double> double2_tuple;
typedef std::map<POSE::atom_type, double2_tuple> atom_type_double2_tuple_map;
typedef boost::tuple<POSE::atom_type, POSE::atom_type> atom_type2_tuple;
typedef std::map<atom_type2_tuple, double2_tuple> atom_type2_tuple_double2_tuple_map;

potential_shared_ptr new_lj_nb_potential();

class lj_nb_potential : public potential_interface {


	friend potential_shared_ptr new_lj_nb_potential();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const{
		return 0;
	}

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const{
		return 0;
	}


	bool init();

	static const double scale_to_joules;// = 4184; //scale kcal to J
	static const double R;// = 8.314472; //J K-1 mol-1
	static const double assumed_temperature;// = 200.0;
	static const double scale_factor;// = (scale_to_joules / (R * assumed_temperature));


protected:
	lj_nb_potential();

	atom_type_double2_tuple_map atomVdwParams;
	atom_type2_tuple_double2_tuple_map atomPairVdwParams;


	void add_lj_params(const POSE::atom_type at1, const double sigma, const double epsilon ){
		atomVdwParams[at1] = double2_tuple(sigma, epsilon);
	}

	void add_lj_params(const POSE::atom_type at1, const POSE::atom_type at2, const double sigma, const double epsilon ){
		atomPairVdwParams[atom_type2_tuple(at1, at2)] = double2_tuple(sigma, epsilon);
		atomPairVdwParams[atom_type2_tuple(at2, at1)] = double2_tuple(sigma, epsilon);
	}

	void calc_pair_vals(){
		for (atom_type_double2_tuple_map::const_iterator it1 = atomVdwParams.begin(); it1 != atomVdwParams.end(); it1++){
			for (atom_type_double2_tuple_map::const_iterator it2 = atomVdwParams.begin(); it2 != atomVdwParams.end(); it2++){

				const double sigma = std::sqrt(it1->second.get<0>() * it2->second.get<0>());
				const double epsilon = std::sqrt(it1->second.get<1>() * it2->second.get<1>());
				add_lj_params(it1->first, it2->first, sigma, epsilon);

			}
		}
	}



	double vdw_atr(const double dist_pow6, const double sigma_pow6, const double epsilon) const{
		const double part_6 = (sigma_pow6 / dist_pow6);

		return  ( - 4.0 * epsilon * part_6);
	}

	double vdw_rep(const double dist_pow12, const double sigma_pow12, const double epsilon) const{
		const double part_12 = ( sigma_pow12 / dist_pow12);

		return ( 4.0 * epsilon * part_12);
	}


	double vdw_atr_grad_mag(const double dist_pow7, const double sigma_pow6, const double epsilon) const{
		const double part_6 = - (sigma_pow6 / dist_pow7);

		return  ( - 6.0 * 4.0 * epsilon *part_6);
	}

	double vdw_rep_grad_mag(const double dist_pow13, const double sigma_pow12, const double epsilon) const{
		const double part_12 = -( sigma_pow12 / dist_pow13);

		return ( 12.0 * 4.0 * epsilon *part_12);
	}

	void get_lj_energy(const PRODART::POSE::META::nb_pair_element& ele, double& atr, double& rep ) const{

		atr = 0;
		rep = 0;

		if (ele.atom1_ptr->isActive() && ele.atom2_ptr->isActive()){
			atom_type2_tuple at_pair(ele.atype1, ele.atype2);
			if (atomPairVdwParams.find(at_pair) != atomPairVdwParams.end()){
				double2_tuple params = atomPairVdwParams.find(atom_type2_tuple(ele.atype1, ele.atype2))->second;
				const double sigma = params.get<0>();
				const double epsilon = params.get<1>();

				const double sigma_pow6 = std::pow(sigma, 6);
				const double dist_pow6 = std::pow(ele.dist_sq, 3);

				atr = vdw_atr(dist_pow6, sigma_pow6, epsilon);//- 4.0 * epsilon *  std::pow(sigma / ele.dist_sq, 3);
				rep = vdw_rep(std::pow(dist_pow6, 2), std::pow(sigma_pow6, 2), epsilon);
			}
		}

	}

	void get_lj_energy_grad(const PRODART::POSE::META::nb_pair_element& ele,
			double& atr,
			double& rep,
			PRODART::UTILS::vector3d_vector& grad,
			const double weight_atr,
			const double weight_rep ) const{

		atr = 0;
		rep = 0;

		if (ele.atom1_ptr->isActive() && ele.atom2_ptr->isActive()){
			atom_type2_tuple at_pair(ele.atype1, ele.atype2);
			if (atomPairVdwParams.find(at_pair) != atomPairVdwParams.end()){

				double2_tuple params = atomPairVdwParams.find(atom_type2_tuple(ele.atype1, ele.atype2))->second;
				const double sigma = params.get<0>();
				const double epsilon = params.get<1>();

				const double sigma_pow6 = std::pow(sigma, 6);
				const double sigma_pow12 = std::pow(sigma_pow6, 2);
				const double dist_pow6 = std::pow(ele.dist_sq, 3);

				atr = vdw_atr(dist_pow6, sigma_pow6, epsilon);//- 4.0 * epsilon *  std::pow(sigma / ele.dist_sq, 3);
				rep = vdw_rep(std::pow(dist_pow6, 2), sigma_pow12, epsilon);


				const PRODART::UTILS::vector3d unit_grad_vec = (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()) / ele.dist;
				const double dist_pow7 = std::pow(ele.dist, 7);
				const double dist_pow13 = std::pow(ele.dist, 13);

				const double grad_mag = (weight_atr * vdw_atr_grad_mag(dist_pow7, sigma_pow6, epsilon))
						+ (weight_rep * vdw_rep_grad_mag(dist_pow13, sigma_pow12, epsilon));

				const PRODART::UTILS::vector3d atm1_grad = grad_mag * unit_grad_vec;
				const PRODART::UTILS::vector3d atm2_grad = -grad_mag * unit_grad_vec;

				const int index1 = ele.atom1_ptr->get_seq_num();
				const int index2 = ele.atom2_ptr->get_seq_num();
				grad[index1] += atm1_grad;
				grad[index2] += atm2_grad;


			}
		}

	}

	void get_lj_energy(const PRODART::POSE::META::nb_ele_vector& vec, double& atr, double& rep ) const{

		atr = 0;
		rep = 0;


		for (PRODART::POSE::META::nb_ele_vector::const_iterator iter = vec.begin(); iter != vec.end(); iter++ ){
			double this_atr = 0, this_rep = 0;
			get_lj_energy(*iter, this_atr, this_rep);
			atr += this_atr;
			rep += this_rep;

		}

		atr *= this->scale_factor;
		rep *= this->scale_factor;

	}

	void get_lj_energy_grad(const PRODART::POSE::META::nb_ele_vector& vec,
			double& atr,
			double& rep,
			PRODART::UTILS::vector3d_vector& grad,
			const double weight_atr,
			const double weight_rep ) const{

		atr = 0;
		rep = 0;


		for (PRODART::POSE::META::nb_ele_vector::const_iterator iter = vec.begin(); iter != vec.end(); iter++ ){
			double this_atr = 0, this_rep = 0;
			get_lj_energy_grad(*iter, this_atr, this_rep, grad, this->scale_factor * weight_atr, this->scale_factor * weight_rep);
			atr += this_atr;
			rep += this_rep;

		}

		atr *= this->scale_factor;
		rep *= this->scale_factor;

	}



};


}
}
}

#endif /* LJ_NB_POTENTIAL_H_ */
