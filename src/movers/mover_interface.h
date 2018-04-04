//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * mover_interface.h
 *
 *  Created on: 1 Mar 2010
 *      Author: jmacdona
 */

#ifndef MOVER_INTERFACE_H_
#define MOVER_INTERFACE_H_

#include "pose/pose.h"
#include "pose_meta/pose_meta_interface.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include <boost/shared_ptr.hpp>
#include "mover_defs.h"
#include <vector>
#include "potentials/potentials_name.h"

namespace PRODART {
namespace POSE {
namespace MOVERS {

typedef std::vector<bool> bool_vector;

class mover_interface;

typedef boost::shared_ptr<mover_interface> mover_shared_ptr;

class mover_interface {


public:

	enum mover_type {
		mt_pure_mover,
		mt_move_set
	};

	virtual bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const = 0;

	void set_rand_num_gen(MTRand::MTRand_shared_ptr ptr){
		rand_gen = ptr;
	}

	virtual void propagate_rand_num_gen(const MTRand::MTRand_shared_ptr ptr){
		rand_gen = ptr;
	}

	//! for composite move sets
	virtual void propagate_start_beta(const double val){
	}
	virtual void propagate_final_beta(const double val){
	}
	virtual void propagate_ca_potential_preset(const std::string str){
	}
	virtual void propagate_bb_potential_preset(const std::string str){
	}
	virtual void propagate_ca_potential_weight(const POTENTIALS::potentials_name&, const double weight){
	}
	virtual void propagate_bb_potential_weight(const POTENTIALS::potentials_name&, const double weight){
	}

	mover_type get_type() const{
		return type;
	}

	void auto_rand_num_assign(){
		if (!rand_gen){
			rand_gen = PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id());
		}
	}

protected:

	MTRand::MTRand_shared_ptr rand_gen;
	mover_type type;

	mover_interface(){
		type = mt_pure_mover;
	}

	virtual ~mover_interface(){
	}

private:




};

}
}
}



#endif /* MOVER_INTERFACE_H_ */
