/*
 * misc_protocols.h
 *
 *  Created on: 18 Nov 2010
 *      Author: jmacdona
 */

#ifndef MISC_PROTOCOLS_H_
#define MISC_PROTOCOLS_H_

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <vector>
#include <iomanip>
#include <fstream>
#include "prodart_env/prodart_env.h"
#include <exception>
#include "pose_meta/ca_pose_meta.h"
#include "pose_meta/bb_pose_meta.h"
#include "pose_meta/frag_classify/fragment_classifier_interface.h"
#include "pose_meta/frag_classify/ca_fragment_classifier.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "utils/line_fit.h"
#include "pose/pose.h"
#include "pose_utils/pose_utils.h"
#include "prodart_env/prodart_env.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "movers/ca_movers/ca_single_dihedral_mover.h"
#include "movers/ca_movers/ca_single_angle_mover.h"
#include "movers/ca_movers/ca_single_crankshaft_mover.h"
#include "pose_meta/ca_pose_meta.h"
#include "potentials/ca_potentials/ca_bump.h"
#include "potentials/bond_potential.h"
#include "potentials/potentials_container.h"
#include "potentials/potentials_factory.h"
#include "simulation/monte_carlo.h"
#include "simulation/ca_mc_protocols/simple_ca_mc_protocol.h"
#include "simulation/minimiser.h"
#include "backbonebuilder/backbone_builder.h"
#include "pose_utils/residue_reconstructor.h"
#include "rotamers/rotamer_db.h"
#include "rotamers/sidechain_builders/sidechain_builder.h"


namespace PRODART {
namespace PROTOCOLS{
namespace MISC{

typedef std::vector<bool> bool_vector;

PRODART::POSE::pose_shared_ptr verified_load_ca_pose( const std::string filename );
PRODART::POSE::pose_shared_ptr verified_load_bb_pose( const std::string filename );
PRODART::POSE::pose_shared_ptr unverified_load_pose( const std::string filename );

PRODART::POSE::three_state_sec_struct_vector read_3state_secs(const std::string filename );

void output_rosetta_match_csts(PRODART::POSE::pose_shared_ptr protein,
		const std::string filename,
		bool use_atom_types = false);
void output_rosetta_match_csts(PRODART::POSE::atom_shared_ptr res1atm3,
		PRODART::POSE::atom_shared_ptr res1atm2,
		PRODART::POSE::atom_shared_ptr res1atm1,
		PRODART::POSE::atom_shared_ptr res2atm1,
		PRODART::POSE::atom_shared_ptr res2atm2,
		PRODART::POSE::atom_shared_ptr res2atm3,
		bool use_atom_types = false);

void output_rosetta_contraints_file_polyALA(PRODART::POSE::pose_shared_ptr protein, std::ostream& output);
void output_rosetta_contraints_file_all(PRODART::POSE::pose_shared_ptr protein, std::ostream& output);
void output_rosetta_contraints_file_all(PRODART::POSE::pose_shared_ptr protein, bool_vector loop_mask, std::ostream& output);
void output_rosetta_sc_sc_contraints_file(PRODART::POSE::pose_shared_ptr protein, bool_vector loop_mask, std::ostream& output);
void trim_terminal_tails(PRODART::POSE::pose_shared_ptr protein, const int max_tail);
void inactivate_unset_GLY_CBs(PRODART::POSE::pose_shared_ptr protein);

void output_rama_histogram_from_raw_phi_psi( std::istream& input,
		std::ostream& output,
		const bool trans = true,
		const bool pro = false);
void output_rama_histogram_from_raw_phi_psi_filelist( std::istream& input,
		std::ostream& output,
		const bool trans = true);

void output_raw_phi_psi_omega(PRODART::POSE::pose_shared_ptr protein, std::ostream& output, std::string label);
void output_raw_phi_psi_omega_from_list(std::istream& input, std::ostream& output);

std::istream& next_ENDMDL(std::istream& input);

PRODART::POSE::pose_shared_ptr get_single_mdl(const int mdl_num, std::istream& input);
void clean_pose(PRODART::POSE::pose_shared_ptr protein, bool remove_bb_h = false);

void to_bb_only(PRODART::POSE::pose_shared_ptr protein);

void graft_motif_residues_to_scaffold(PRODART::POSE::pose_shared_ptr scaffold, PRODART::POSE::pose_shared_ptr motif);

// make a pose with single idealised ALA residue chain A and CA at (0,0,0), CB at (1.5461,0,0), C in the XY plane, etc...
PRODART::POSE::pose_shared_ptr make_centred_ALA_residue();
PRODART::POSE::pose_shared_ptr make_centred_trans_pept_bond();
PRODART::POSE::pose_shared_ptr make_centred_cis_pept_bond();
PRODART::POSE::pose_shared_ptr copy_backbone_segment(PRODART::POSE::const_pose_shared_ptr protein,
		const int start,
		const int end,
		const bool renumber = true,
		const bool minimal_bb_set = true,
		const bool include_sc = false);
bool output_centred_backbone_segments(PRODART::POSE::pose_shared_ptr protein,
		const int fraglen,
		std::ostream& output,
		const bool not_gly = true,
		const bool not_pro = true);
bool output_centred_backbone_segments(std::istream& filelist,
		const int fraglen,
		std::ostream& output,
		const bool not_gly = true,
		const bool not_pro = true);

bool output_quality_filtered_list(std::istream& filelist,
		std::ostream& output);

//! best to avoid
PRODART::POSE::pose_shared_ptr get_single_chain_pose_byChainID(PRODART::POSE::const_pose_shared_ptr pose_, char chainID);
//! prefered function
PRODART::POSE::pose_shared_ptr get_single_chain_pose_byChainNum(PRODART::POSE::const_pose_shared_ptr pose_, const int chainnum);

bool process_pdb_file_for_rosetta_match(std::string& pose_file,
		const std::string& output_dir,
		const int min_len,
		const int max_len);
bool process_list_for_rosetta_match(std::istream& filelist,
		const std::string& output_dir,
		const int min_len,
		const int max_len);


bool output_design_quality(std::istream& filelist,
		 std::ostream& output);

bool output_seq_diff( std::ostream& output,
		PRODART::POSE::pose_shared_ptr ref,
		PRODART::POSE::pose_shared_ptr protein);

bool output_seq3( std::ostream& output,
		PRODART::POSE::pose_shared_ptr protein);

bool output_seq1( std::ostream& output,
		PRODART::POSE::pose_shared_ptr protein);

class except_bool_bad_parse: public std::exception
{
  virtual const char* what() const throw()
  {
    return "ERROR: couldn't parse string into bool_vector";
  }
};
inline std::vector<bool> make_bool_vec(std::string str){
	std::vector<bool> ret_vec(str.size(), false);

	for ( unsigned int i = 0 ; i < str.size(); i++){
		if (str[i] == '1'){
			ret_vec[i] = true;
		}
		else if (str[i] == '0'){
			ret_vec[i] = false;
		}
		else {
			throw except_bool_bad_parse();
		}
	}
	return ret_vec;

}

boost::tuple<double, double, double> get_max_alpha_beta_loop_N_CA_C_movement(const PRODART::POSE::pose_shared_ptr ref,
		const PRODART::POSE::pose_shared_ptr protein);

//! returns number of cis-prolines in scaffold turned into non-proline residues
int get_cisPro_to_nonPro_count(const PRODART::POSE::pose_shared_ptr ref,
		const PRODART::POSE::pose_shared_ptr protein);



PRODART::POSE::atom_shared_ptr_vector get_atom_subset(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& residue_loop_mask, PRODART::POSE::atom_type_vector &vec);

PRODART::POSE::atom_shared_ptr_vector get_atom_subset(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& residue_loop_mask, std::vector<PRODART::POSE::BBAtomType> &vec);

void output_raw_bond_lengths_angles(PRODART::POSE::const_pose_shared_ptr protein, std::ostream& output);

void output_raw_bond_lengths_angles( std::istream& input_filelist,
		std::ostream& output);

//! returns the maximum backbone atom deviation in angstroms from the motif
double get_max_motif_atom_deviation(const int start_res, PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		double& rtn_rmsd);


void scan_motif_atom_deviation(POSE::double_vector& min_vals,
		std::vector<std::string>& mot_res,
		PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const bool use_rmsd = false);

//! outputs the lowest maximum backbone atom deviation in angstroms from the motif by sliding motif along the protein
void scan_motif_atom_deviation(std::ostream& output,
		PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const bool use_rmsd = false);



void scan_motif_atom_deviation(std::ostream& output,
		 std::istream& input_filelist,
		PRODART::POSE::const_pose_shared_ptr motif,
		const bool use_rmsd = false);




//! outputs the lowest maximum backbone atom deviation in angstroms from the motif by sliding motif along the protein
void scan_motif_atom_deviation_annot_seq(std::ostream& output,
		PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd = false);


void scan_motif_atom_deviation_annot_seq(std::ostream& output,
		 std::istream& input_filelist,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd = false);

//! outputs all fragments
void scan_motif_atom_deviation_output_frags(PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd = false);


void scan_motif_atom_deviation_output_frags(std::istream& input_filelist,
		PRODART::POSE::const_pose_shared_ptr motif,
		const double cutoff,
		const bool use_rmsd = false);


//! inverts chain CAs then rebuilds backbone
PRODART::POSE::pose_shared_ptr invert_chain_direction(PRODART::POSE::const_pose_shared_ptr protein);

std::map< int,std::map< int, boost::tuple<double, double, double> > > make_psicov_map(PRODART::POSE::const_pose_shared_ptr protein,
		const double cb_cutoff,
		const int min_seq_sep = 5);
std::map< int,std::map< int, boost::tuple<double, double, double> > > parse_psicov_file(std::istream& psicov_input);

std::ostream& output_psicov_contact_map(PRODART::POSE::const_pose_shared_ptr protein,
		const double cb_cutoff,
		const int min_seq_sep ,
		std::ostream& output);


//seq sep based on residue numbers only
std::ostream& compare_to_psicov_contact_map(PRODART::POSE::const_pose_shared_ptr protein,
		const double cb_cutoff,
		const double prob_cutoff,
		const int min_seq_sep,
		std::istream& psicov_input,
		std::ostream& output);

std::ostream& output_contactmap_rstfile(PRODART::POSE::pose_shared_ptr protein,
		const double cutoff,
		std::ostream& output);

PRODART::POSE::pose_shared_ptr make_extended_ca_chain(const PRODART::POSE::residue_type_vector seq,
		const char chainID,
		const int start_res,
		const double random_component = 0.1);

void extend_ca_chain(PRODART::POSE::pose_shared_ptr  protein);

PRODART::POSE::pose_shared_ptr make_dummy_ca_pose(const PRODART::POSE::residue_type_vector seq,
		const char chainID,
		const int start_res = 1,
		const double rand_coord = 0);

bool thread_sequence(const PRODART::POSE::residue_type_vector seq,
		PRODART::POSE::pose_shared_ptr pose_);


std::vector<boost::tuple< std::string, char, std::string, char  > > parse_residue_mapping(std::istream& input_residue_mapping);

std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> get_const_atom_mapping(
		const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::const_pose_shared_ptr protein_to_move,
		std::vector<boost::tuple< std::string, char, std::string, char  > > residue_map,
		bool verbose = true);

//! calculate RMSD given an arbitrary atom mapping - can be used with any atoms - then translates and rotates the supplied pose (should but doesn't have to correspond to 'second' in the atom map) - NOT rigorously tested yet
double get_rmsd_superpose_atmmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::pose_shared_ptr protein_to_move,
		std::istream& input_atom_mapping);


//! calculate RMSD given an arbitrary residue mapping - can be used with any atoms - then translates and rotates the supplied pose (should but doesn't have to correspond to 'second' in the atom map) - NOT rigorously tested yet
double get_rmsd_superpose_resmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::pose_shared_ptr protein_to_move,
		std::istream& input_residue_mapping,
		bool verbose = true);

//! calculate RMSD given an arbitrary residue mapping without any superposition - can be used with any atoms - NOT rigorously tested yet
double get_rmsd_no_superpose_resmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		const PRODART::POSE::const_pose_shared_ptr protein_to_move,
		std::istream& input_residue_mapping,
		bool verbose = true);

//! calculate vector of RMSDs given an arbitrary residue mapping without any superposition - can be used with any atoms - NOT rigorously tested yet
std::vector<double> get_rmsd_no_superpose_trajectory_resmap(const PRODART::POSE::const_pose_shared_ptr ref_protein,
		std::istream& traj_input,
		std::istream& input_residue_mapping,
		bool verbose = false);


//! output Rosetta format all-atom distance restraints
std::ostream& output_rosetta_distance_restraints_all_atom(PRODART::POSE::const_pose_shared_ptr protein,
		const double max_dist_cutoff,
		const double stddev,
		std::ostream& output);

//! output Rosetta format all-atom harmonic coordinate restraints
std::ostream& output_rosetta_coord_restraints_all_atom(PRODART::POSE::const_pose_shared_ptr protein,
		const double stddev,
		std::ostream& output);

bool search_replace_atom_types(const PRODART::POSE::pose_shared_ptr pose_, std::istream& search_replace_mapping);


}
}
}




#endif /* MISC_PROTOCOLS_H_ */
