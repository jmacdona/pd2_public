/*
 * ca_single_angle_mover.h
 *
 *  Created on: 9 Sep 2010
 *      Author: jmacdona
 */

#ifndef CA_SINGLE_ANGLE_MOVER_H_
#define CA_SINGLE_ANGLE_MOVER_H_



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


class ca_single_angle_uni_dist_mover;

typedef boost::shared_ptr<ca_single_angle_uni_dist_mover> ca_single_angle_uni_dist_mover_shared_ptr;

//! uniform distribution mover
class ca_single_angle_uni_dist_mover : public mover_interface {

	friend ca_single_angle_uni_dist_mover_shared_ptr new_ca_single_angle_uni_dist_mover(const bool is_forwards,
			const double max_angle,
			const int ang_num);

public:

	bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const;

private:

	ca_single_angle_uni_dist_mover();
	ca_single_angle_uni_dist_mover(const bool is_forwards,
			const double max_angle,
			const int ang_num);
	void init();

	double max_rot_ang;
	int angle_num;
	bool isFwds;


};



ca_single_angle_uni_dist_mover_shared_ptr new_ca_single_angle_uni_dist_mover(const bool is_forwards,
		const double max_angle,
		const int ang_num);

move_set_shared_ptr ca_single_angle_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_angle,
		const PRODART::POSE::MOVERS::bool_vector& allowed_residues);
move_set_shared_ptr ca_single_angle_uni_dist_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const double max_angle);






}
}
}
}











#endif /* CA_SINGLE_ANGLE_MOVER_H_ */
