/*
 * sse_axes.h
 *
 *  Created on: Mar 12, 2012
 *      Author: jmacdon
 */

#ifndef SSE_AXES_H_
#define SSE_AXES_H_

#include "pose/pose.h"
#include "utils/vector3d.h"
#include "utils/line_fit.h"
#include "qcprot.h"
#include "pose_basic_kine.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>




namespace PRODART {
namespace POSE_UTILS {


void calculate_helix_axis(const PRODART::UTILS::vector3d& c1,
		const PRODART::UTILS::vector3d& c2,
		const PRODART::UTILS::vector3d& c3,
		const PRODART::UTILS::vector3d& c4,
		PRODART::UTILS::vector3d& start,
		PRODART::UTILS::vector3d& end);

//! full extent helix calculation
void calculate_helix_axis_long(const PRODART::UTILS::vector3d& c1,
		const PRODART::UTILS::vector3d& c2,
		const PRODART::UTILS::vector3d& c3,
		const PRODART::UTILS::vector3d& c4,
		PRODART::UTILS::vector3d& start,
		PRODART::UTILS::vector3d& end);

void calculate_strand_axis(const PRODART::UTILS::vector3d& c1,
		const PRODART::UTILS::vector3d& c2,
		const PRODART::UTILS::vector3d& c3,
		const PRODART::UTILS::vector3d& c4,
		PRODART::UTILS::vector3d& start,
		PRODART::UTILS::vector3d& end);

//! full extent strand calculation
void calculate_strand_axis_long(const PRODART::UTILS::vector3d& c1,
		const PRODART::UTILS::vector3d& c2,
		const PRODART::UTILS::vector3d& c3,
		const PRODART::UTILS::vector3d& c4,
		PRODART::UTILS::vector3d& start,
		PRODART::UTILS::vector3d& end);

//! full extent helix calculation
bool get_helix_end_points_long(const PRODART::UTILS::vector3d_vector& ca_coords,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec);

//! full extent strand calculation
bool get_strand_end_points_long(const PRODART::UTILS::vector3d_vector& ca_coords,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec);

//! full extent helix calculation
bool get_helix_end_points_long(const PRODART::POSE::const_pose_shared_ptr protein,
		const int start_residue_index,
		const int end_residue_index,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec);

//! full extent strand calculation
bool get_strand_end_points_long(const PRODART::POSE::const_pose_shared_ptr protein,
		const int start_residue_index,
		const int end_residue_index,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec);

}
}
#endif /* SSE_AXES_H_ */
