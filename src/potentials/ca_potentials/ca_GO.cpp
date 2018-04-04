/*
 * ca_GO.cpp
 *
 *  Created on: 3 Oct 2013
 *      Author: jmacdona
 */
#include "ca_GO.h"
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
namespace CA{




potential_shared_ptr new_ca_GO(){
	potential_shared_ptr ptr(new ca_GO());
	return ptr;
}

ca_GO::ca_GO() : upper_limit(12), min_seq_sep(5) {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_GO"));
}






bool ca_GO::init(){

	return true;
}


double ca_GO::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();

	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(pose_meta_);

	PRODART::POSE::META::nb_ele_vector& pair_list = ca_meta_dat->get_ca_pair_list();

	PRODART::POSE::META::nb_ele_vector::iterator iter;

	double  total_energy = 0;

	//int invalid_count = 0;

	const double penalty = 1.0;

	if (pose_meta_->are_go_restraints_set() ){

		for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
			if (iter->dist < upper_limit ){

				//const int res1 = std::min(iter->res_num1, iter->res_num2);
				//const int res2 = std::max(iter->res_num1, iter->res_num2);
				if (abs(iter->res_num2 - iter->res_num1) >= min_seq_sep){
					upper_lower_bonded_pair_element ele = ca_meta_dat->get_ca_GO_restraint(iter->res_num1, iter->res_num2);
					if (ele.is_valid){
						if (iter->dist < ele.equilibrium_dist_upper && iter->dist > ele.equilibrium_dist_lower){
							//cout << "VALID" << endl;

							//assume ele.half_bond_const is positive so:
							const double inc = -ele.half_bond_const;
							total_energy += inc;
							/*
					cout << "GOOD contact: res1:" << iter->res_num1 << " res2:"
							<< iter->res_num2 << " inc:" << inc
							<< " dist:" << iter->dist
							<< " ele.equilibrium_dist_lower:" <<ele.equilibrium_dist_lower
							<< " ele.equilibrium_dist_upper:" <<ele.equilibrium_dist_upper
							<< endl;
							 */
						}

					}
					else {
						//Penalise non-native contact
						//assume ele.half_bond_const is positive so penalty:
						const double inc = penalty;//ele.half_bond_const;
						total_energy += inc;
					}
				}
			}
		}


	}

	return energies_map.add_energy_component(name_vector[0], total_energy);

}

double ca_GO::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	return get_energy(pose_meta_, energies_map);//energies_map.add_energy_component(name_vector[0], total_energy);

}















}
}
}
}





