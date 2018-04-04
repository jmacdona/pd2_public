/*
 * bb_14lj.cpp
 *
 *  Created on: 22 Nov 2010
 *      Author: jmacdona
 */

#include "bb_14_15lj.h"

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

potential_shared_ptr new_bb_14_15lj(){
	potential_shared_ptr ptr(new bb_14_15lj());
	return ptr;
}

bb_14_15lj::bb_14_15lj(){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_14_lj_atr"));
	this->name_vector.push_back(potentials_name("bb_14_lj_rep"));
	this->name_vector.push_back(potentials_name("bb_15_lj_atr"));
	this->name_vector.push_back(potentials_name("bb_15_lj_rep"));

	this->add_lj_params(atom_type("C"), 3.750, 0.105);
	this->add_lj_params(atom_type("O"), 2.960 , 0.210 );
	this->add_lj_params(atom_type("N"), 3.250, 0.170);
	this->add_lj_params(atom_type("CA"), 3.800 , 0.080);
	this->add_lj_params(atom_type("CB"), 3.910, 0.160);

	this->calc_pair_vals();



	/*
	atomVdwParams[C].diameter = 3.750;
	atomVdwParams[C].epsilon = 0.105;

	atomVdwParams[O].diameter = 2.960;
	atomVdwParams[O].epsilon = 0.210;

	atomVdwParams[N].diameter = 3.250;
	atomVdwParams[N].epsilon = 0.170;

	atomVdwParams[H].diameter = 0.000;
	atomVdwParams[H].epsilon = 0.000;

	atomVdwParams[CA].diameter = 3.800;
	atomVdwParams[CA].epsilon = 0.080;

	atomVdwParams[CB].diameter = 3.910;
	atomVdwParams[CB].epsilon = 0.160;
	*/

}

double bb_14_15lj::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();


	double en14atr = 0, en14rep = 0;
	double en15atr = 0, en15rep = 0;

	const PRODART::POSE::META::nb_ele_vector& lj14 = bb_meta_dat->get_bonded_14_list();
	this->get_lj_energy(lj14,en14atr, en14rep );
	const PRODART::POSE::META::nb_ele_vector& lj15 = bb_meta_dat->get_bonded_15_list();
	this->get_lj_energy(lj15,en15atr, en15rep);

	return energies_map.add_energy_component(name_vector[0], en14atr)
			+ energies_map.add_energy_component(name_vector[1], en14rep)
			+ energies_map.add_energy_component(name_vector[2], en15atr)
			+ energies_map.add_energy_component(name_vector[3], en15rep);

}

double bb_14_15lj::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight14atr = energies_map.get_weight(this->name_vector[0]);
	const double weight14rep = energies_map.get_weight(this->name_vector[1]);
	const double weight15atr = energies_map.get_weight(this->name_vector[2]);
	const double weight15rep = energies_map.get_weight(this->name_vector[3]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();

	double en14atr = 0, en14rep = 0;
	double en15atr = 0, en15rep = 0;

	const PRODART::POSE::META::nb_ele_vector& lj14 = bb_meta_dat->get_bonded_14_list();
	this->get_lj_energy_grad(lj14,en14atr, en14rep, grad, weight14atr, weight14rep  );
	const PRODART::POSE::META::nb_ele_vector& lj15 = bb_meta_dat->get_bonded_15_list();
	this->get_lj_energy_grad(lj15,en15atr, en15rep, grad, weight15atr, weight15rep );

	//assert(is_vector3d_vector_valid(grad));


	return energies_map.add_energy_component(name_vector[0], en14atr)
			+ energies_map.add_energy_component(name_vector[1], en14rep)
			+ energies_map.add_energy_component(name_vector[2], en15atr)
			+ energies_map.add_energy_component(name_vector[3], en15rep);
}


bool bb_14_15lj::init(){

	return true;
}




}
}
}
}

