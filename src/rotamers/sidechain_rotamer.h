/*
 * sidechain_rotamer.h
 *
 *  Created on: 14 Feb 2011
 *      Author: jmacdona
 */

#ifndef SIDECHAIN_ROTAMER_H_
#define SIDECHAIN_ROTAMER_H_

#include "pose/pose.h"
#include <map>
#include <boost/tuple/tuple.hpp>
#include "pose_utils/pose_basic_kine.h"
#include "rotamer_entry.h"

namespace PRODART {

namespace POSE_UTILS{
class residue_reconstructor;
}

namespace ROTAMERS {

class sidechain_rotamer;

typedef boost::shared_ptr<sidechain_rotamer> sidechain_rotamer_shared_ptr;



class sidechain_rotamer : public PRODART::POSE::sidechain{

	friend sidechain_rotamer_shared_ptr new_sidechain_rotamer();
	friend class PRODART::POSE_UTILS::residue_reconstructor;

private:



protected:
	sidechain_rotamer();
	sidechain_rotamer(const sidechain_rotamer&);


	std::map<POSE::BBAtomType, bool> backbone_mask;
	std::map<POSE::BBAtomType, boost::tuple<UTILS::vector3d, POSE::atom_type> > backbone_changes;
	PRODART::POSE::residue_type res_type;

	const POSE::prot_backbone_map*  bb_map;

public:

	~sidechain_rotamer(){
	}

	void apply() const;

	void set_residue_type(PRODART::POSE::residue_type ty){
		res_type = ty;
	}

	void set_backbone_mask(std::map<POSE::BBAtomType, bool> &mask){
		backbone_mask = mask;
	}

	void set_bacbone_changes(std::map<POSE::BBAtomType, boost::tuple<UTILS::vector3d, POSE::atom_type> >& changes){
		backbone_changes = changes;
	}

	bool set_chi(const unsigned int chi_index, const double rad){
		return POSE_UTILS::KINE::set_chi_angle_forwards(this->_this.lock(),chi_index, rad, true );
	}
	bool set_chi_inverse(const unsigned int chi_index, const double rad){
		return POSE_UTILS::KINE::set_chi_angle_backwards(this->_this.lock(),chi_index, rad, true );
	}



};

sidechain_rotamer_shared_ptr new_sidechain_rotamer();

}
}

#endif /* SIDECHAIN_ROTAMER_H_ */
