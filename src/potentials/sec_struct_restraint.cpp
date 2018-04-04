/*
 * sec_struct_restraint.cpp
 *
 *  Created on: 27 Jan 2011
 *      Author: jmacdona
 */

#include "sec_struct_restraint.h"

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

potential_shared_ptr new_sec_struct_restraint(){
	potential_shared_ptr ptr(new sec_struct_restraint());
	return ptr;
}

sec_struct_restraint::sec_struct_restraint() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("sec_struct_restraint"));
}



bool sec_struct_restraint::init(){

	return true;
}

double sec_struct_restraint::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const PRODART::POSE::four_state_sec_struct_vector& conf_class = pose_meta_->get_conf_class();
	const PRODART::POSE::three_state_sec_struct_vector& sec_struct = pose_meta_->get_sec_struct();
	const PRODART::POSE::four_state_sec_struct_vector& rsts = pose_meta_->get_sec_struct_restraint_list();
	const double_vector& weights = pose_meta_->get_sec_struct_restraint_weight_list();

	double total_energy = 0;

	for (unsigned int i = 0; i < rsts.size(); i++){
		const four_state_sec_struct ss_rst = rsts[i];
		if (ss_rst != ss4_UNDEF){
			const four_state_sec_struct cclass = conf_class[i];
			const three_state_sec_struct secs = sec_struct[i];
			const double weight = weights[i];

			if (ss_rst == cclass){
				total_energy += -0.5 * weight;
			}

			if (ss_rst == ss4_HELIX && secs == ss3_HELIX){
				total_energy += -0.5 * weight;
			}
			else if (ss_rst == ss4_STRAND && secs == ss3_STRAND){
				total_energy += -0.5 * weight;
			}
			else if (ss_rst == ss4_OTHER && secs == ss3_OTHER){
				total_energy += -0.5 * weight;
			}


		}
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double sec_struct_restraint::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	return this->get_energy(pose_meta_, energies_map);
}


}
}
}
