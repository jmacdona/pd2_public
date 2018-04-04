//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose_utils.h
 *
 *  Created on: Feb 19, 2010
 *      Author: jmacdona
 */

#ifndef POSE_UTILS_H_
#define POSE_UTILS_H_

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
#include <boost/tuple/tuple.hpp>




namespace PRODART {
namespace POSE_UTILS {

/*
 *
 * NOTE that this has been declared in pose/residue.h :
 *
enum three_state_sec_struct{
	secs_OTHER = 0,
	secs_HELIX = 1,
	secs_STRAND = 2
*
*
*/


//! calculate line segment for the given residue range and fitting to CA atoms; startVec and endVec return the end points
bool get_line_CA_segment(const PRODART::POSE::const_pose_shared_ptr protein,
		const int start_residue_index,
		const int end_residue_index,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec);



//! calculate carbon-alpha 1-to-1 RMSD - NOT rigorously tested yet
double get_ca_rmsd(const PRODART::POSE::const_pose_shared_ptr protein1,
		const PRODART::POSE::const_pose_shared_ptr protein2);

//! calculate carbon-alpha 1-to-1 RMSD and superimpose second pose on first pose - NOT rigorously tested yet but seems to work
double get_ca_rmsd_superpose(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::pose_shared_ptr protein_to_move);

//! calculate carbon-alpha RMSD using alignment given by residue_mapping - NOT rigorously tested yet
double get_ca_rmsd(const PRODART::POSE::const_pose_shared_ptr protein1,
		const PRODART::POSE::const_pose_shared_ptr protein2,
		const std::map<int,int> residue_mapping);

//! calculate carbon-alpha RMSD using alignment given by residue_mapping and superpose second pose on ref pose - NOT rigorously tested yet
double get_ca_rmsd_superpose(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::pose_shared_ptr protein_to_move,
		const std::map<int,int> residue_mapping);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - NOT rigorously tested yet
double get_rmsd(const std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - NOT rigorously tested yet
double get_msd(const std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - NOT rigorously tested yet
double get_msd(const std::vector< boost::tuple <PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr, double > > atom_atom_wt_mapping,
		const PRODART::UTILS::vector3d& CoM1, const PRODART::UTILS::vector3d& CoM2);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - NOT rigorously tested yet
double get_rmsd(const std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - NOT rigorously tested yet
double get_msd(const std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - then translates and rotates the supplied pose (should but doesn't have to correspond to 'second' in the atom map) - NOT rigorously tested yet
double get_rmsd_superpose(const PRODART::POSE::pose_shared_ptr protein_to_move,
		const std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - then translates and rotates the supplied pose (should but doesn't have to correspond to 'second' in the atom map) - NOT rigorously tested yet
double get_rmsd_superpose(const PRODART::POSE::pose_shared_ptr protein_to_move,
		const std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - then translates and rotates the supplied pose (should but doesn't have to correspond to 'second' in the atom map) - NOT rigorously tested yet
double get_rmsd_superpose(const PRODART::POSE::pose_shared_ptr protein_to_move,
		const std::vector<boost::tuple<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr, double> >  atom_atom_wt_mapping,
		const PRODART::UTILS::vector3d& CoM_ref, const PRODART::UTILS::vector3d& CoM_to_move);

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - then translates and rotates the supplied atoms (should but doesn't have to correspond to 'second' in the atom map) - NOT rigorously tested yet
double get_rmsd_superpose(const PRODART::POSE::atom_shared_ptr_vector atoms_to_move,
		const std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping);












double get_fixed_position_rmsd(PRODART::POSE::atom_shared_ptr_vector& vec1, PRODART::POSE::atom_shared_ptr_vector& vec2);

double get_fixed_position_rmsd(const std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> &atom_mapping);

double get_fixed_position_rmsd_auto_loop(const PRODART::POSE::const_pose_shared_ptr pose1,
		const PRODART::POSE::pose_shared_ptr pose2,
		const double tol = 0.001);

//! centre of mass based in CA atoms
PRODART::UTILS::vector3d get_ca_CoM(const PRODART::POSE::const_pose_shared_ptr protein);



//! get radius of gyration given a residue range
double get_ca_radgyr(const PRODART::POSE::const_pose_shared_ptr protein,
		const int startRes,
		const int endRes);

//! get radius of gyration given a chain index
double get_ca_radgyr(const PRODART::POSE::const_pose_shared_ptr protein,
		const int chain_num);

//! quick and dirty addition of backbone hydrogen to a residue
UTILS::vector3d get_ideal_HN(const PRODART::POSE::pose_shared_ptr protein,
		const int residue_number,
		bool &is_defined);

//! quick and dirty addition of backbone hydrogen to a residue
bool quick_add_HN(const PRODART::POSE::pose_shared_ptr protein,
		const int residue_number,
		const bool overwrite);


//! quick and dirty addition of backbone hydrogen to a whole pose
bool quick_add_HN(const PRODART::POSE::pose_shared_ptr protein,
		const bool overwrite);

UTILS::vector3d get_rough_ideal_CB(const PRODART::POSE::const_pose_shared_ptr protein,
		const int residue_number);

//! quick and dirty addition of CB to a residue
bool quick_add_CB(const PRODART::POSE::pose_shared_ptr protein,
		const int residue_number,
		const bool overwrite = false);

//! quick and dirty addition of  CB to a whole pose including GLYs
bool quick_add_CB(const PRODART::POSE::pose_shared_ptr protein,
		const bool overwrite = false);

//! quick and dirty addition of  CB to a whole pose not including GLYs
bool quick_add_CB_notGLY(const PRODART::POSE::pose_shared_ptr protein,
		const bool overwrite = false);

//! find lowest and highest bounding coordinates values
void get_ca_bounding_box(const PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::UTILS::vector3d& lowest, PRODART::UTILS::vector3d& highest);


inline POSE::four_state_sec_struct get_phi_psi_omega_sector(const double Phi, const double Psi, const double omega){

	if ( omega < (UTILS::PI / 2.0) && omega > -(UTILS::PI / 2.0)) {
		//is cis
		return POSE::ss4_CIS;
	}
	else if (Phi > -UTILS::PI && Phi < -((20.0 * UTILS::PI)/180.0) &&
                    Psi > -((90.0 * UTILS::PI)/180.0) && Psi < -((10.0 * UTILS::PI)/180.0) ){
            //alpha helix

            return POSE::ss4_HELIX;
    }
    else if ( Phi >= -UTILS::PI && Phi < -((20.0 * UTILS::PI)/180.0) && (
                    (Psi <= UTILS::PI && Psi > ((20.0 * UTILS::PI)/180.0))
                    || (Psi >= -UTILS::PI && Psi < -((170.0 * UTILS::PI)/180.0)) )
    ){
            //beta strand

            return POSE::ss4_STRAND;
    }
    //other
    return POSE::ss4_OTHER;
}

POSE::four_state_sec_struct_vector get_phi_psi_omega_sector(POSE::const_pose_shared_ptr pose_);

void to_homopolymer(const POSE::pose_shared_ptr pose_, const POSE::residue_type rt);

bool validate_bb_residue(const PRODART::POSE::const_residue_shared_ptr res);

bool validate_bb_pose(const PRODART::POSE::const_pose_shared_ptr protein,
		const double bond_tol = 0.1);

bool validate_ca_pose(const PRODART::POSE::const_pose_shared_ptr protein);

//! read single chain FASTA - ignores lines beginning with ">"
PRODART::POSE::residue_type_vector read_fasta_single_chain(const std::string filename);

bool strip_bb_hydrogens(const POSE::pose_shared_ptr pose_);

bool search_replace_atom_types(const POSE::pose_shared_ptr pose_, const std::map<POSE::atom_type, POSE::atom_type> s_r_map);


}
}








#endif /* POSE_UTILS_H_ */
