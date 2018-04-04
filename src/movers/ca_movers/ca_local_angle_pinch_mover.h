/*
 * ca_local_angle_pinch_mover.h
 *
 *  Created on: 28 Sep 2010
 *      Author: jmacdona
 */

#ifndef CA_LOCAL_ANGLE_PINCH_MOVER_H_
#define CA_LOCAL_ANGLE_PINCH_MOVER_H_

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


class ca_local_angle_pinch_mover;

typedef boost::shared_ptr<ca_local_angle_pinch_mover> ca_local_angle_pinch_mover_shared_ptr;

//! uniform distribution mover
class ca_local_angle_pinch_mover : public mover_interface {

	friend ca_local_angle_pinch_mover_shared_ptr new_ca_local_angle_pinch_mover(const double max_angle,
			const int bondNum);

public:

	bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const;

private:

	ca_local_angle_pinch_mover();
	ca_local_angle_pinch_mover(const double max_angle,
			const int bondNum);
	void init();

	double max_rot_ang;
	int angleNum;


};



ca_local_angle_pinch_mover_shared_ptr new_ca_local_angle_pinch_mover(const double max_angle,
			const int bondNum);

move_set_shared_ptr ca_local_angle_pinch_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_dist,
		const PRODART::POSE::MOVERS::bool_vector& allowed_residues);

move_set_shared_ptr ca_local_angle_pinch_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_dist);






}
}
}
}





#endif /* CA_LOCAL_ANGLE_PINCH_MOVER_H_ */
