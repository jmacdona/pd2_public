//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * ca_single_dihedral_mover.h
 *
 *  Created on: 2 Mar 2010
 *      Author: jmacdona
 */

#ifndef CA_SINGLE_DIHEDRAL_MOVER_H_
#define CA_SINGLE_DIHEDRAL_MOVER_H_

#include "movers/mover_interface.h"
#include "movers/move_set.h"



#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include <boost/shared_ptr.hpp>

namespace PRODART {
namespace POSE {
namespace MOVERS {
namespace CA {

class ca_single_dihedral_uni_dist_mover;

typedef boost::shared_ptr<ca_single_dihedral_uni_dist_mover> ca_single_dihedral_uni_dist_mover_shared_ptr;

//! uniform distribution mover
class ca_single_dihedral_uni_dist_mover : public mover_interface {

	friend ca_single_dihedral_uni_dist_mover_shared_ptr new_ca_single_dihedral_uni_dist_mover(const bool is_forwards,
			const double max_angle,
			const int ang_num,
			const int largest_move_resnum);

public:

	bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const;

private:

	ca_single_dihedral_uni_dist_mover();
	ca_single_dihedral_uni_dist_mover(const bool is_forwards,
			const double max_angle,
			const int ang_num,
			const int largest_move_resnum);
	void init();

	bool isFwds;
	double max_rot_ang;
	int angle_num;
	int big_move_resnum;



};

ca_single_dihedral_uni_dist_mover_shared_ptr new_ca_single_dihedral_uni_dist_mover(const bool is_forwards,
		const double max_angle,
		const int ang_num,
		const int largest_move_resnum);
move_set_shared_ptr ca_single_dihedral_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data, const double max_angle,
		const PRODART::POSE::MOVERS::bool_vector& allowed_residues);
move_set_shared_ptr ca_single_dihedral_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data, const double max_angle);

}
}
}
}


#endif /* CA_SINGLE_DIHEDRAL_MOVER_H_ */
