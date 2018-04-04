/*
 * ca_pseudo_cb_pos.cpp
 *
 *  Created on: 22 Feb 2011
 *      Author: jmacdona
 */

#ifndef CA_PSEUDO_CB_POS_CPP_
#define CA_PSEUDO_CB_POS_CPP_


#include "ca_pseudo_cb_pos.h"

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


ca_pseudo_cb_pos_shared_ptr new_ca_pseudo_cb_pos(){
	ca_pseudo_cb_pos_shared_ptr ptr(new ca_pseudo_cb_pos());
	return ptr;
}

ca_pseudo_cb_pos::ca_pseudo_cb_pos() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_pseudo_cb_pos"));

	b_is_disposable = true;
}






double ca_pseudo_cb_pos::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	double total_score = 0;

	const double ideal_CA_CB_bond = 1.526;

	for (std::map<int, UTILS::vector3d>::const_iterator it = resnum_ca_pos_map.begin(); it != resnum_ca_pos_map.end(); it++){
		const int resnum = it->first;


		if (((resnum-1) >= 0) && ((resnum + 1) < _pose->get_residue_count())){

			const vector3d cb_pos = it->second;

			const vector3d ca_m1 = _pose->get_bb_atom(POSE::CA, resnum-1)->get_coords();
			const vector3d ca = _pose->get_bb_atom(POSE::CA, resnum)->get_coords();
			const vector3d ca_p1 = _pose->get_bb_atom(POSE::CA, resnum+1)->get_coords();

			const vector3d v1 = ca - ca_m1;
			const vector3d v2 = ca - ca_p1;
			const vector3d unit = (v1 + v2) / (v1 + v2).mod();

			const vector3d pseu_cb = ca + (ideal_CA_CB_bond * unit);

			total_score += (cb_pos - pseu_cb).mod();
		}

	}


	return energies_map.add_energy_component(this->name_vector[0], total_score);
}

//TODO
double ca_pseudo_cb_pos::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{

	return get_energy(_pose_meta, energies_map);

}


void ca_pseudo_cb_pos::add_rst(const int resnum, const UTILS::vector3d vec){
	resnum_ca_pos_map[resnum] = vec;
}

bool ca_pseudo_cb_pos::init(){

	return true;
}

//! does it provide energies labeled with this name
/*
bool ca_pseudo_cb_pos::provides(const potentials_name& query_name) const{
	if (query_name == name){
		return true;
	}

	return false;
}
*/










}
}
}
}


#endif /* CA_PSEUDO_CB_POS_CPP_ */
