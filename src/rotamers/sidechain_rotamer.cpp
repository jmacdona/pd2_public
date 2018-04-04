/*
 * sidechain_rotamer.cpp
 *
 *  Created on: 14 Feb 2011
 *      Author: jmacdona
 */

#include "sidechain_rotamer.h"
#include "pose_utils/residue_reconstructor.h"

using namespace PRODART::POSE;
using namespace PRODART;
using namespace UTILS;

namespace PRODART {
namespace ROTAMERS {

//const prot_backbone_map*  sidechain_rotamer::bb_map = prot_backbone_map::Instance();



void sidechain_rotamer::apply() const{
	residue_shared_ptr res_to_apply = res.lock();
	assert(res_to_apply);
	pose_shared_ptr pose_ = res_to_apply->get_pose();
	assert(pose_);
	assert(_this.lock());
	res_to_apply->swap_sidechain(_this.lock());

	{
		std::map<POSE::BBAtomType, bool>::const_iterator it;
		for (it = backbone_mask.begin(); it != backbone_mask.end(); it++){
			atom_shared_ptr at = res_to_apply->get_bb_atom(it->first);
			at->setSet(it->second);
			at->setActive(it->second);
		}
	}
	{
		std::map<POSE::BBAtomType, boost::tuple<UTILS::vector3d, POSE::atom_type> >::const_iterator it;
		for (it = backbone_changes.begin(); it != backbone_changes.end(); it++){
			atom_shared_ptr at = res_to_apply->get_bb_atom(it->first);
			at->setActive(true);
			at->setSet(true);
			at->set_coords(it->second.get<0>());
			at->set_type(it->second.get<1>());
		}
	}
	res_to_apply->set_type(res_type);
	_this.lock()->update_chi_atom_cache();
	pose_->index();
}


sidechain_rotamer::sidechain_rotamer() : sidechain(), bb_map(prot_backbone_map::Instance()){
	const std::vector<BBAtomType> full_bb_map = bb_map->get_full_bb_atom_list();
	backbone_mask.clear();
	for (unsigned int i = 0; i < full_bb_map.size(); i++){
		backbone_mask[full_bb_map[i]] = false;
	}

}

sidechain_rotamer::sidechain_rotamer(const sidechain_rotamer& old_sc) : sidechain(old_sc), bb_map(prot_backbone_map::Instance()){

}

sidechain_rotamer_shared_ptr new_sidechain_rotamer(){
	sidechain_rotamer_shared_ptr ptr(new sidechain_rotamer);
	ptr->_this = ptr;
	return ptr;
}



}
}

