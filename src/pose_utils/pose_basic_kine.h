//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose_basic_kine.h
 *
 *  Created on: 23 Feb 2010
 *      Author: jmacdona
 */

#ifndef POSE_BASIC_KINE_H_
#define POSE_BASIC_KINE_H_

#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
//#include "qcprot.h"
//#include "pose_utils.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>
#include <cassert>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>




namespace PRODART {
namespace POSE_UTILS {
namespace KINE {



//! translate all atoms by translation vector
void translate_fa(const PRODART::POSE::pose_shared_ptr protein,
		const PRODART::UTILS::vector3d& trans);

//! translate supplied atom vector by translation vector
void translate_fa(const PRODART::POSE::atom_shared_ptr_vector atoms_to_move,
		const PRODART::UTILS::vector3d& trans);

//! applies rotation matrix to all atoms - assumes the pose has already been centered
void apply_rotation_matrix_fa(const PRODART::POSE::pose_shared_ptr protein,
		const PRODART::UTILS::rot_matrix& rot);

//! applies rotation matrix to supplied atom vector - assumes the pose has already been centered
void apply_rotation_matrix_fa(const PRODART::POSE::atom_shared_ptr_vector atoms_to_move,
		const PRODART::UTILS::rot_matrix& rot);

void rotate_dihedral_forwards_ca(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives = false);

void rotate_dihedral_backwards_ca(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives = false);

double rotate_dihedral_forwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives = false);

double rotate_dihedral_backwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives = false);

double crankshaft_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int resNum1,
		const int resNum2,
		const double angle_rot,
		const bool run_through_inactives = false);


double rotate_angle_forwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives = false);

double rotate_angle_backwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives = false);

bool set_chi_angle_forwards(const PRODART::POSE::atom_shared_ptr_vector& chi_def,
		const PRODART::POSE::atom_shared_ptr_vector& chi_fwd,
		const double chi,
		const bool run_through_inactives = false);

bool set_chi_angle_backwards(const PRODART::POSE::atom_shared_ptr_vector& chi_def,
		const PRODART::POSE::atom_shared_ptr_vector& chi_bwd,
		const double chi,
		const bool run_through_inactives = false);

bool set_chi_angle_forwards(PRODART::POSE::sidechain_shared_ptr sc,
		const unsigned int chi_index,
		const double chi,
		const bool run_through_inactives = false);

bool set_chi_angle_backwards(PRODART::POSE::sidechain_shared_ptr sc,
		const unsigned int chi_index,
		const double chi,
		const bool run_through_inactives = false);

}
}
}

#endif /* POSE_BASIC_KINE_H_ */
