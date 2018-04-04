//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pdb_utils.cpp
 *
 *  Created on: 5 Mar 2010
 *      Author: jmacdona
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include "utils/line_fit.h"
#include "pose/residue_type.h"
#include "pose/atom_type.h"
#include "pose/atom.h"
#include "pose/residue.h"
#include "pose/pose.h"
#include "pose_utils/pose_utils.h"
#include "prodart_env/prodart_env.h"
#include "pose_meta/ca_pose_meta.h"
#include "protocols/protocols.h"


using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::pow;
using std::sqrt;
using std::asin;
using std::string;

using boost::lexical_cast;
using boost::split;
using boost::is_any_of;


typedef std::vector< std::string >  split_vector_type;


using PRODART::UTILS::PI;
using PRODART::POSE::atom_type;
using PRODART::UTILS::vector3d;

using namespace PRODART::ENV;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::META;
using namespace PRODART::PROTOCOLS;
using namespace std;

//! very simple pdb utility functions
int main( int argc, char *argv[] ) {


	//options_manager &opt_mgr = prodart_env::Instance()->get_options_manager();
	//PRODART::ENV::list_options_full_format(cout);
	PRODART::ENV::register_option("help",bool(false),"print options");
	PRODART::ENV::register_option("renumber",bool(false),"renumber residues");
	PRODART::ENV::register_option("renumber:start",int(1),"renumber residues start resnum");
	PRODART::ENV::register_option("renumber:cyclic",bool(false),"cyclic renumber residues");
	PRODART::ENV::register_option("renumber:cyclic:period",int(5),"cyclic renumber period");
	PRODART::ENV::register_option("get_rmsd",bool(false),"calculate RMSD");
	PRODART::ENV::register_option("get_rmsd_no_superpose",bool(false),"calculate RMSD");
	PRODART::ENV::register_option("strip_bb_hydrogens",bool(false),"strip all backbone hydrogen atoms");
	PRODART::ENV::register_option("renumber_chainIDs",bool(false),"renumber chain IDs");
	PRODART::ENV::register_option("fill_missing_chainIDs",bool(false),"fill blank chain IDs with an unused character");


	PRODART::ENV::register_option("get_rmsd_auto_loop",bool(false),"calculate RMSD");
	PRODART::ENV::register_option("get_rmsd_auto_loop:tol",double(0.001),"tolerance");
	PRODART::ENV::register_option("ca2main:get_rmsd_special",bool(false),"get RMSD for ca2main method - special for Ben's ca2main paper. Ignore first NCONCON and last CONCONCO");
	PRODART::ENV::register_option("ca2main:get_rmsd_sc",bool(false),"get RMSD for ca2main method - special for Ben's ca2main paper. Only sidechains not CBs");
	PRODART::ENV::register_option("get_ca_secs",bool(false),"calculate sec struct from CA coords using fragments and hbonding patterns");
	PRODART::ENV::register_option("get_secs",bool(false),"calculate sec struct from BB coords using phi/psi and hbonding patterns");
	PRODART::ENV::register_option("get_hot_spots",bool(false),"print backbone hot spots");
	PRODART::ENV::register_option("get_phi_psi_omega",bool(false),"print phi psi angles");
	PRODART::ENV::register_option("get_single_mdl",bool(false),"get single model from PDB trajectory file");
	PRODART::ENV::register_option("clean_pdb",bool(false),"clean PDB occupancies etc");
	PRODART::ENV::register_option("get_backbone",bool(false),"get backbone atoms only in PDB");

	PRODART::ENV::register_option("get_ca_only",bool(false),"get CA PDB");
	PRODART::ENV::register_option("graft_motif_residues_to_scaffold",bool(false),"graft motif residues to nearest scaffold residues");
	PRODART::ENV::register_option("pdb_utils:motif_pdb", string(""),"input file");
	PRODART::ENV::register_option("add_CBs",bool(false),"add carbon betas to PDB");
	PRODART::ENV::register_option("add_CBs:overwrite_existing",bool(false),"overwrite existing CBs");

	PRODART::ENV::register_option("mdl_num", int(1), "model number in PDB trajectory - starts from 1");
	PRODART::ENV::register_option("make_rosetta_constraints_file",bool(false),"output_rosetta_contraints_file_polyALA");
	PRODART::ENV::register_option("make_rosetta_constraints_file_all",bool(false),"output_rosetta_contraints_file_all");
	PRODART::ENV::register_option("make_rosetta_sc_sc_constraints_file",bool(false),"output_rosetta_sc_sc_contraints_file");
	PRODART::ENV::register_option("make_rosetta_all_atom_constraints_file",bool(false),"Outputs all atom distance constraints file for Rosetta using the output_rosetta_distance_restraints_all_atom function");
	PRODART::ENV::register_option("make_rosetta_all_atom_coord_cst_file",bool(false),"Outputs all atom coordinate constraints file for Rosetta using the output_rosetta_distance_restraints_all_atom function");
	PRODART::ENV::register_option("make_rosetta_match_cst_file",bool(false),"output_rosetta_match_csts");
	PRODART::ENV::register_option("make_rosetta_match_cst_file:rm_atoms", string(""),"output_rosetta_match_csts");
	PRODART::ENV::register_option("make_phi_psi_hist",bool(false),"output_phi_psi_hist");
	PRODART::ENV::register_option("make_phi_psi_hist:input",string(""),"input file");
	PRODART::ENV::register_option("make_phi_psi_hist:list",string(""),"input file list");
	PRODART::ENV::register_option("make_phi_psi_hist:trans",bool(true),"trans conformation, if not cis");
	PRODART::ENV::register_option("make_phi_psi_hist:pro",bool(false),"allow prolines as well");
	PRODART::ENV::register_option("trim_terminal_tails",bool(false),"trim_terminal_tails");
	PRODART::ENV::register_option("make_loop_mask",bool(false),"outputs residue loop mask");
	PRODART::ENV::register_option("make_loop_mask:resnum",string(""),"residue_number");
	PRODART::ENV::register_option("chainid",char('A'),"chain ID");

	PRODART::ENV::register_option("process_pdbs_for_rosetta_match",bool(false),"process_pdbs_for_rosetta_match");
	PRODART::ENV::register_option("process_pdbs_for_rosetta_match:min_len",int(1),"process_pdbs_for_rosetta_match");
	PRODART::ENV::register_option("process_pdbs_for_rosetta_match:max_len",int(numeric_limits<int>::max()),"process_pdbs_for_rosetta_match");
	PRODART::ENV::register_option("output_dir",string("./"),"output directory");

	PRODART::ENV::register_option("output_design_quality",bool(false),"output design quality scores");
	PRODART::ENV::register_option("output_design_quality:input",string(""),"input list. 2 columns - 1) design pdb 2) original scaffold");
	PRODART::ENV::register_option("output_design_quality:output",string(""),"output file");



	PRODART::ENV::unhide_option("loop_model:loop_mask");

	PRODART::ENV::register_option("make_dummy_pdb",bool(false),"dummy PDB");
	PRODART::ENV::register_option("make_dummy_pdb:seq",string("A"),"1 letter seq");
	PRODART::ENV::register_option("make_dummy_pdb:start_resnum",int(1),"start resnum");
	PRODART::ENV::register_option("make_dummy_pdb:dummy_chain",char('A'),"start resnum");
	PRODART::ENV::register_option("make_dummy_pdb:rand",double(1),"rand");

	PRODART::ENV::register_option("residue_mapping",string(""),"residue mapping file");
	PRODART::ENV::register_option("atom_mapping",string(""),"atom mapping file");

	PRODART::ENV::register_option("make_centered_bb_frags",bool(false),"makes backbone fragments centred on central CA atom at (0,0,0), CB at (1.5461,0,0), C in the XY plane, etc");
	PRODART::ENV::register_option("make_centered_bb_frags:bb_frag_len", int(5), "fragment length - should probably use odd numbers");
	PRODART::ENV::register_option("make_centered_bb_frags:no_gly",bool(true),"no GLY residues in central position");
	PRODART::ENV::register_option("make_centered_bb_frags:no_pro",bool(true),"no PRO residues in central position");

	PRODART::ENV::register_option("output_filtered_list",bool(false),"filter input pdb filelist for quality");
	PRODART::ENV::register_option("output_seq_diff",bool(false),"output sequence diff");

	PRODART::ENV::register_option("get_seq3",bool(false),"get  3 letter sequence list");
	PRODART::ENV::register_option("get_seq",bool(false),"get  3 letter sequence list");


	PRODART::ENV::register_option("get_raw_angles_bonds",bool(false),"print bond angles and bond lengths");

	PRODART::ENV::register_option("scan_repeat_motif",bool(false),"output min max deviation values");
	PRODART::ENV::register_option("scan_repeat_motif:use_rmsd",bool(false),"use backbone RMSD instead of atom deviation");

	PRODART::ENV::register_option("scan_repeat_motif:output_annot_seq",bool(false),"output min max deviation values");
	PRODART::ENV::register_option("scan_repeat_motif:cutoff",double(1.0),"use backbone RMSD instead of atom deviation");
	PRODART::ENV::register_option("scan_repeat_motif:output_frags",bool(false),"output pdb fragments");



	PRODART::ENV::register_option("get_chain",bool(false),"output chain");

	PRODART::ENV::register_option("invert_chain_direction",bool(false),"invert chain direction");

	PRODART::ENV::register_option("get_psicov_contact_map",bool(false),"outputs psicov format contact map");
	PRODART::ENV::register_option("compare_with_psicov_contact_map",bool(false),"outputs psicov format contact map");
	PRODART::ENV::register_option("psicov_map",string(""),"input psicov contact map");
	PRODART::ENV::register_option("psicov_cutoff",double(0.5),"psicov prob cutoff");
	PRODART::ENV::register_option("psicov_min_seq_sep",int(5),"psicov prob cutoff");

	PRODART::ENV::register_option("make_contact_map_rstfile",bool(false),"make Go contact map rstfile");
	PRODART::ENV::register_option("make_contact_map_rstfile:cutoff",double(12.0),"distance cutoff");


	PRODART::ENV::register_option("thread_sequence",bool(false),"thread sequence onto PDB");
	PRODART::ENV::register_option("fasta_in",string(""), "input FASTA format sequence file");

	PRODART::ENV::register_option("pose:io:pdb:no_check_prev_atomid",bool(false),"disable check atom id is incremented");


	PRODART::ENV::register_option("search_replace_atom_names",bool(false),"search and replace atom types");
	PRODART::ENV::register_option("search_replace_atom_names:mapping", string(""),"search and replace atom type mapping file");


	if (!prodart_env::Instance()->init(argc, argv)){
		cerr << "initialisation errors" << endl;
		return -1;
	}

	if (PRODART::ENV::get_option_value<bool>("help") == true){
		PRODART::ENV::list_options_full_format(cout);
		return 0;
	}

	if (PRODART::ENV::is_set("")){
		cerr << "ERROR: command line free options are not supported" << endl;
		return -1;
	}

	const bool input_pdb_is_set = PRODART::ENV::is_set("pose:io:pdb:i");
	const bool output_pdb_is_set = PRODART::ENV::is_set("pose:io:pdb:o");
	const bool ref_pdb_is_set = PRODART::ENV::is_set("pose:io:pdb:r");


	PRODART::POSE::pose_shared_ptr in_pdb = PRODART::POSE::new_pose();
	ifstream protein_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:i").c_str(), ios::in);
	if (input_pdb_is_set){
		if (protein_file.is_open()){
			in_pdb->loadPdb(protein_file, !PRODART::ENV::get_option_value<bool>("pose:io:pdb:no_check_prev_atomid"));
			in_pdb->set_label(PRODART::ENV::get_option_value<string>("pose:io:pdb:i"));
		}
		else {
			cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:i")
			                                           << endl;
			return -1;
		}
		protein_file.close();
	}
	/*
	else{
		cerr << "ERROR: no pdb input file set\n"
			 << endl;
		return -1;
	}
	*/

	PRODART::POSE::pose_shared_ptr ref_pdb = PRODART::POSE::new_pose();
	ifstream ref_protein_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:r").c_str(), ios::in);
	if (ref_pdb_is_set){
		if (ref_protein_file.is_open()){
			ref_pdb->loadPdb(ref_protein_file, !PRODART::ENV::get_option_value<bool>("pose:io:pdb:no_check_prev_atomid"));
			ref_pdb->set_label(PRODART::ENV::get_option_value<string>("pose:io:pdb:r"));
		}
		else {
			cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:r")
			                                           << endl;
			return -1;
		}
		ref_protein_file.close();
	}

	ofstream output_pdb_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:o").c_str(), ios::out);
	if (output_pdb_is_set){
		if (!output_pdb_file.is_open()){
			cerr << "ERROR: can't open output file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:o")
			                                           << endl;
			return -1;
		}
	}

	if (input_pdb_is_set){
		if (PRODART::ENV::get_option_value<bool>("renumber") == true){
			//! renumber residues
			if (PRODART::ENV::get_option_value<bool>("renumber:cyclic") == false){
				in_pdb->renumber_residues(ENV::get_option_value<int>("renumber:start"));
			}
			else {
				in_pdb->cyclic_renumber_residues(ENV::get_option_value<int>("renumber:start"),
						1,
						ENV::get_option_value<int>("renumber:cyclic:period"),
						ENV::get_option_value<int>("renumber:cyclic:period"));
			}
			if (output_pdb_is_set){
				in_pdb->outputPdb(output_pdb_file);
			}
			else {
				in_pdb->outputPdb(cout);
			}
			return 0;
		}
		else if (PRODART::ENV::get_option_value<bool>("renumber_chainIDs") == true){
			in_pdb->renumber_chainIDs();
			if (output_pdb_is_set){
				in_pdb->outputPdb(output_pdb_file);
			}
			else {
				in_pdb->outputPdb(cout);
			}
			return 0;
		}
		else if (PRODART::ENV::get_option_value<bool>("fill_missing_chainIDs") == true){
			in_pdb->auto_assign_missing_chainIDs();
			if (output_pdb_is_set){
				in_pdb->outputPdb(output_pdb_file);
			}
			else {
				in_pdb->outputPdb(cout);
			}
			return 0;
		}
		else if (PRODART::ENV::get_option_value<bool>("get_chain") == true){
			pose_shared_ptr pose_ = PROTOCOLS::MISC::get_single_chain_pose_byChainID(in_pdb, ENV::get_option_value<char>("chainid"));
			if (output_pdb_is_set){
				pose_->outputPdb(output_pdb_file);
			}
			else {
				pose_->outputPdb(cout);
			}
			return 0;

		}
		else if (PRODART::ENV::get_option_value<bool>("search_replace_atom_names") == true
				&& PRODART::ENV::is_set("search_replace_atom_names:mapping")
				&& output_pdb_is_set){
			ifstream at_map(PRODART::ENV::get_option_value<string>("search_replace_atom_names:mapping").c_str(), ios::in);
			if (at_map.is_open()){
				PROTOCOLS::MISC::search_replace_atom_types(in_pdb, at_map);
				in_pdb->outputPdb(output_pdb_file);
			}
			else {
				cerr << "ERROR: couldn't open mapping file" << endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("get_rmsd") == true
				&& PRODART::ENV::is_set("residue_mapping")
				&& !PRODART::ENV::is_set("atom_mapping")
				&& output_pdb_is_set
				&& ref_pdb_is_set){
			ifstream res_map(PRODART::ENV::get_option_value<string>("residue_mapping").c_str(), ios::in);
			if (res_map.is_open()){
				cout << PRODART::PROTOCOLS::MISC::get_rmsd_superpose_resmap(ref_pdb,
						in_pdb,
						res_map) << endl;
				in_pdb->outputPdb(output_pdb_file);




			}
			else {
				cerr << "ERROR: couldn't open residue_mapping file" << endl;
			}
		}
		else if (PRODART::ENV::get_option_value<bool>("get_rmsd_no_superpose") == true
				&& PRODART::ENV::is_set("residue_mapping")
				&& !PRODART::ENV::is_set("atom_mapping")
				&& ref_pdb_is_set){
			ifstream res_map(PRODART::ENV::get_option_value<string>("residue_mapping").c_str(), ios::in);
			if (res_map.is_open()){

				cout << PRODART::PROTOCOLS::MISC::get_rmsd_no_superpose_resmap(ref_pdb,
						in_pdb,
						res_map) << endl;




			}
			else {
				cerr << "ERROR: couldn't open residue_mapping file" << endl;
			}
		}
		else if (PRODART::ENV::get_option_value<bool>("ca2main:get_rmsd_special") == true
				&& ref_pdb_is_set){


				PRODART::PROTOCOLS::CA2MAIN::ca2main_get_rmsd_special(ref_pdb, in_pdb);


		}
		else if (PRODART::ENV::get_option_value<bool>("ca2main:get_rmsd_sc") == true
				&& ref_pdb_is_set){


				PRODART::PROTOCOLS::CA2MAIN::ca2main_get_rmsd_sidechain(ref_pdb, in_pdb);


		}
		else if (PRODART::ENV::get_option_value<bool>("get_rmsd") == true
				&& !PRODART::ENV::is_set("residue_mapping")
				&& PRODART::ENV::is_set("atom_mapping")
				&& output_pdb_is_set
				&& ref_pdb_is_set){
			ifstream atm_map(PRODART::ENV::get_option_value<string>("atom_mapping").c_str(), ios::in);
			if (atm_map.is_open()){
				cout << PRODART::PROTOCOLS::MISC::get_rmsd_superpose_atmmap(ref_pdb,
						in_pdb,
						atm_map) << endl;
				in_pdb->outputPdb(output_pdb_file);




			}
			else {
				cerr << "ERROR: couldn't open residue_mapping file" << endl;
			}
		}
		else if (PRODART::ENV::get_option_value<bool>("get_rmsd_auto_loop") == true
				&& input_pdb_is_set
				&& ref_pdb_is_set){

			cout << PRODART::POSE_UTILS::get_fixed_position_rmsd_auto_loop(ref_pdb,
					in_pdb,
					PRODART::ENV::get_option_value<double>("get_rmsd_auto_loop:tol")) << endl;

		}
		else if (PRODART::ENV::get_option_value<bool>("get_ca_secs") == true){
			ca_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_ca_pose_meta(in_pdb);
			cout << pose_meta_->get_sec_struct() << endl;
		}
		else if (PRODART::ENV::get_option_value<bool>("get_secs") == true){
			POSE_UTILS::quick_add_HN(in_pdb,false);
			for (int i = 0; i < in_pdb->get_residue_count(); i++){
				atom_shared_ptr at = in_pdb->get_bb_atom(POSE::H, i);
				at->set_type(atom_type("H"));
			}
			bb_pose_meta_shared_ptr pose_meta_= PRODART::POSE::META::new_bb_pose_meta(in_pdb);
			cout << pose_meta_->get_sec_struct() << endl;
		}
		else if (PRODART::ENV::get_option_value<bool>("get_hot_spots") == true){
			POSE_UTILS::quick_add_HN(in_pdb,false);
			LOOP::print_bb_hot_spots(in_pdb, cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_rosetta_constraints_file") == true){
			MISC::output_rosetta_contraints_file_polyALA(in_pdb, cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_rosetta_all_atom_constraints_file") == true){
			MISC::output_rosetta_distance_restraints_all_atom(in_pdb, ENV::get_option_value<double>("rosetta_contraints:cst_all_atom_cutoff"), 1, cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_rosetta_all_atom_coord_cst_file") == true){
			// TODO
			//asdfsadfasdfasdfasdf;
			//rosetta_contraints:cst_all_atom_cutoff
			MISC::output_rosetta_coord_restraints_all_atom(in_pdb, 1, cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_rosetta_constraints_file_all") == true && !ENV::is_set("loop_model:loop_mask")){
			MISC::output_rosetta_contraints_file_all(in_pdb, cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_rosetta_constraints_file_all") == true && ENV::is_set("loop_model:loop_mask")){
			LOOP::bool_vector loop_mask = PROTOCOLS::MISC::make_bool_vec(PRODART::ENV::get_option_value<string>("loop_model:loop_mask"));
			MISC::output_rosetta_contraints_file_all(in_pdb, loop_mask, cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_rosetta_sc_sc_constraints_file") == true && ENV::is_set("loop_model:loop_mask")){
			LOOP::bool_vector loop_mask = PROTOCOLS::MISC::make_bool_vec(PRODART::ENV::get_option_value<string>("loop_model:loop_mask"));
			MISC::output_rosetta_sc_sc_contraints_file(in_pdb, loop_mask, cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_rosetta_match_cst_file") == true){
			// *************************************************************
			if (ENV::is_set("make_rosetta_match_cst_file:rm_atoms")){
				MISC::output_rosetta_match_csts(in_pdb, ENV::get_option_value<string>("make_rosetta_match_cst_file:rm_atoms"), true);
			}
			else {
				cerr << "ERROR: need the option 'make_rosetta_match_cst_file:rm_atoms' to be set" << endl;
			}
		}
		else if (ENV::get_option_value<bool>("make_contact_map_rstfile") == true){
			// TODO
			MISC::output_contactmap_rstfile(in_pdb, ENV::get_option_value<double>("make_contact_map_rstfile:cutoff"), cout);
		}
		else if (PRODART::ENV::get_option_value<bool>("trim_terminal_tails") == true
				&& output_pdb_is_set){
			MISC::trim_terminal_tails(in_pdb,2);
			in_pdb->outputPdb(output_pdb_file);
		}
		else if (PRODART::ENV::get_option_value<bool>("make_loop_mask") == true){
			residue_shared_ptr rp = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("make_loop_mask:resnum"), PRODART::ENV::get_option_value<char>("chainid"));
			CA2MAIN::bool_vector mask = CA2MAIN::make_residue_loop_mask(in_pdb, rp->get_internal_residue_index(), rp->get_internal_residue_index());
			for (unsigned int i = 0; i < mask.size(); i++){
				cout << mask[i];
			}
			cout << endl;
		}
		else if (PRODART::ENV::get_option_value<bool>("clean_pdb") == true){
			if (output_pdb_is_set){
				MISC::clean_pose(in_pdb, true);
				if (PRODART::ENV::get_option_value<bool>("get_backbone") == true){
					MISC::to_bb_only(in_pdb);
				}
				in_pdb->outputPdb(output_pdb_file);
			}
			else{
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("thread_sequence") == true){
			cout << "thread_sequence: reading files" << endl;
			if (output_pdb_is_set && ENV::is_set("fasta_in")){
				PRODART::POSE::residue_type_vector seqvec = POSE_UTILS::read_fasta_single_chain(ENV::get_option_value<string>("fasta_in"));
				MISC::thread_sequence(seqvec, in_pdb);
				in_pdb->outputPdb(output_pdb_file);
			}
			else{
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("invert_chain_direction") == true){
			if (output_pdb_is_set){
				pose_shared_ptr npose = MISC::invert_chain_direction(in_pdb);
				npose->outputPdb(output_pdb_file);
			}
			else{
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (ENV::get_option_value<bool>("get_psicov_contact_map") == true
				&& output_pdb_is_set){
			MISC::output_psicov_contact_map(in_pdb,
					8,
					ENV::get_option_value<int>("psicov_min_seq_sep"),
					output_pdb_file);
		}
		else if (ENV::get_option_value<bool>("compare_with_psicov_contact_map") == true
				&& ENV::is_set("psicov_map")){
			ifstream psicov_input(PRODART::ENV::get_option_value<string>("psicov_map").c_str(), ios::in);

			if (psicov_input.is_open()){
				MISC::compare_to_psicov_contact_map(in_pdb,
						8,
						ENV::get_option_value<double>("psicov_cutoff"),
						ENV::get_option_value<int>("psicov_min_seq_sep"),
						psicov_input,
						std::cout);
			}
			else {
				cerr << "ERROR: could not open psicov_map file\n";
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("get_backbone") == true){
			if (output_pdb_is_set){
				MISC::to_bb_only(in_pdb);
				if (PRODART::ENV::get_option_value<bool>("clean_pdb") == true){
					MISC::clean_pose(in_pdb, true);
				}
				in_pdb->outputPdb(output_pdb_file);

			}
			else{
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("add_CBs") == true){
			if (output_pdb_is_set){
				POSE_UTILS::quick_add_CB_notGLY(in_pdb, PRODART::ENV::get_option_value<bool>("add_CBs:overwrite_existing"));
				in_pdb->outputPdb(output_pdb_file);
			}
			else{
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("get_ca_only") == true){
			if (output_pdb_is_set){
				in_pdb->outputPdb(output_pdb_file, true);
			}
			else{
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("get_phi_psi_omega") == true
				&& PRODART::ENV::is_set("pose:io:pdb:list:l")){

			if (output_pdb_is_set){
				ifstream ppo_in(PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l").c_str(), ios::in);
				if (!ppo_in.is_open()){
					cerr << "ERROR: can't open input file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l")
				                                        						   << endl;
					return -1;
				}
				else {
					cout << "calculating phi psi omega dihedrals..." << endl;
					MISC::output_raw_phi_psi_omega_from_list(ppo_in, output_pdb_file);
				}
			}
			else {
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("get_phi_psi_omega") == true){

					if (output_pdb_is_set){
							cout << "calculating phi psi omega dihedrals..." << endl;
							MISC::output_raw_phi_psi_omega(in_pdb, output_pdb_file,
									PRODART::ENV::get_option_value<string>("pose:io:pdb:i"));
					}
					else {
						cout << "calculating phi psi omega dihedrals..." << endl;
						MISC::output_raw_phi_psi_omega(in_pdb, cout,
								PRODART::ENV::get_option_value<string>("pose:io:pdb:i"));
					}

				}
		else if (PRODART::ENV::get_option_value<bool>("get_raw_angles_bonds") == true
				&& PRODART::ENV::is_set("pose:io:pdb:list:l")){
			if (output_pdb_is_set){
				ifstream ppo_in(PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l").c_str(), ios::in);
				if (!ppo_in.is_open()){
					cerr << "ERROR: can't open intput file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l")
				                                        						   << endl;
					return -1;
				}
				else {
					cout << "calculating bond angles and bond lengths..." << endl;
					MISC::output_raw_bond_lengths_angles(ppo_in, output_pdb_file);
				}
			}
			else {
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("scan_repeat_motif") == true
				&& PRODART::ENV::is_set("pose:io:pdb:list:l")){
			if (ENV::is_set("pdb_utils:motif_pdb") && output_pdb_is_set){

				PRODART::POSE::pose_shared_ptr motif_pdb = PRODART::POSE::new_pose();
				ifstream mot_protein_file(PRODART::ENV::get_option_value<string>("pdb_utils:motif_pdb").c_str(), ios::in);

				if (mot_protein_file.is_open()){
					motif_pdb->loadPdb(mot_protein_file);
					motif_pdb->set_label(PRODART::ENV::get_option_value<string>("pdb_utils:motif_pdb"));
				}
				else {
					cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pdb_utils:motif_pdb")
						                                        		   << endl;
					return -1;
				}
				mot_protein_file.close();

				ifstream ppo_in(PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l").c_str(), ios::in);
				if (!ppo_in.is_open()){
					cerr << "ERROR: can't open intput file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l")
				                                        						   << endl;
					return -1;
				}
				else {

					cout << "loaded motif with " << motif_pdb->get_chain_count() << " chains." << endl;
					for (int i = 0 ; i < motif_pdb->get_chain_count(); i++ ){
						cout << motif_pdb->get_chain(i)->getChainID();
					}
					cout << endl;

					cout << "scanning for motif..." << endl;
					if (!ENV::get_option_value<bool>("scan_repeat_motif:output_annot_seq")
							&& !ENV::get_option_value<bool>("scan_repeat_motif:output_frags")){
						MISC::scan_motif_atom_deviation(output_pdb_file, ppo_in, motif_pdb,
								PRODART::ENV::get_option_value<bool>("scan_repeat_motif:use_rmsd"));
					}
					else if (ENV::get_option_value<bool>("scan_repeat_motif:output_annot_seq")
							&& !ENV::get_option_value<bool>("scan_repeat_motif:output_frags")) {
						MISC::scan_motif_atom_deviation_annot_seq(output_pdb_file, ppo_in, motif_pdb,
								PRODART::ENV::get_option_value<double>("scan_repeat_motif:cutoff"),
								PRODART::ENV::get_option_value<bool>("scan_repeat_motif:use_rmsd"));
					}
					else if(!ENV::get_option_value<bool>("scan_repeat_motif:output_annot_seq")
							&& ENV::get_option_value<bool>("scan_repeat_motif:output_frags")){
						MISC::scan_motif_atom_deviation_output_frags( ppo_in, motif_pdb,
								PRODART::ENV::get_option_value<double>("scan_repeat_motif:cutoff"),
								PRODART::ENV::get_option_value<bool>("scan_repeat_motif:use_rmsd"));
					}
					else {
						cerr << "ERROR: bad combinations of options" << endl;
					}
				}
			}
			else {
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if (PRODART::ENV::get_option_value<bool>("graft_motif_residues_to_scaffold") == true){
			if ( ENV::is_set("pdb_utils:motif_pdb") && output_pdb_is_set){
				PRODART::POSE::pose_shared_ptr motif_pdb = PRODART::POSE::new_pose();
				ifstream mot_protein_file(PRODART::ENV::get_option_value<string>("pdb_utils:motif_pdb").c_str(), ios::in);

				if (mot_protein_file.is_open()){
					motif_pdb->loadPdb(mot_protein_file);
					motif_pdb->set_label(PRODART::ENV::get_option_value<string>("pdb_utils:motif_pdb"));
				}
				else {
					cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pdb_utils:motif_pdb")
						                                        		   << endl;
					return -1;
				}
				mot_protein_file.close();

				MISC::graft_motif_residues_to_scaffold(in_pdb, motif_pdb );
				in_pdb->outputPdb(output_pdb_file);

			}
			else {
				cerr << "ERROR: not enough options set"
						<< endl;
			}

		}
		else if ( PRODART::ENV::get_option_value<bool>("output_seq_diff") == true
				&& ref_pdb_is_set){
			MISC::output_seq_diff(cout, ref_pdb, in_pdb);
		}
		else if (PRODART::ENV::get_option_value<bool>("get_seq3") == true){
			MISC::output_seq3(cout, in_pdb);
		}
		else if (PRODART::ENV::get_option_value<bool>("get_seq") == true){
			MISC::output_seq1(cout, in_pdb);
		}
		else {
			cerr << "Nothing to do it's up to you" << endl;
		}
	}
	else {
		/*
		 * ***************************************************
		 * *******************NO input PDB options************
		 * ***************************************************
		 */

		if(PRODART::ENV::get_option_value<bool>("get_single_mdl") == true
				&& output_pdb_is_set && PRODART::ENV::is_set("pose:io:pdb:traj:t") ){

			ifstream traj_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:traj:t").c_str(), ios::in);

			if (traj_file.is_open()){
				pose_shared_ptr pose_ = MISC::get_single_mdl( PRODART::ENV::get_option_value<int>("mdl_num"), traj_file);
				if (pose_->get_residue_count() > 0){
					MISC::clean_pose(pose_);
					if (ENV::get_option_value<bool>("strip_bb_hydrogens") == true){
						POSE_UTILS::strip_bb_hydrogens(pose_);
					}
					pose_->outputPdb(output_pdb_file);
				}

			}
			else {
				cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:traj:t")
				                                        						   << endl;
				return -1;
			}



		}
		else if (PRODART::ENV::get_option_value<bool>("get_rmsd_no_superpose") == true
				&& PRODART::ENV::is_set("pose:io:pdb:traj:t")
				&& PRODART::ENV::is_set("residue_mapping")
				&& !PRODART::ENV::is_set("atom_mapping")
				&& ref_pdb_is_set){
			ifstream res_map(PRODART::ENV::get_option_value<string>("residue_mapping").c_str(), ios::in);
			ifstream traj_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:traj:t").c_str(), ios::in);

			if (traj_file.is_open()){
				if (res_map.is_open()){

					std::vector<double>  rmsds = PRODART::PROTOCOLS::MISC::get_rmsd_no_superpose_trajectory_resmap(ref_pdb,
							traj_file,
							res_map,
							false);

					/*
					for (int ii = 0; ii < rmsds.size(); ii++){
						cout << (ii+1) << "\t" << rmsds[ii] << endl;
					}
					*/



				}
				else {
					cerr << "ERROR: couldn't open residue_mapping file" << endl;
				}
			}
			else {
				cerr << "ERROR: couldn't open trajectory file" << endl;
			}
		}
		else if (PRODART::ENV::get_option_value<bool>("make_phi_psi_hist") == true){
			if (PRODART::ENV::is_set("make_phi_psi_hist:input")
					&& output_pdb_is_set
					&& !PRODART::ENV::is_set("make_phi_psi_hist:list") ){
				ifstream ppo_in(PRODART::ENV::get_option_value<string>("make_phi_psi_hist:input").c_str(),
						ios::in);

				if (!ppo_in.is_open()){
					cerr << "ERROR: can't open intput file: " << PRODART::ENV::get_option_value<string>("make_phi_psi_hist:input")
					                                        				   << endl;
					return -1;
				}
				else {
					cout << "making histogram" << endl;
					MISC::output_rama_histogram_from_raw_phi_psi(ppo_in,
							output_pdb_file,
							PRODART::ENV::get_option_value<bool>("make_phi_psi_hist:trans"));
				}

			}
			else if (!PRODART::ENV::is_set("make_phi_psi_hist:input")
					&& output_pdb_is_set
					&& PRODART::ENV::is_set("make_phi_psi_hist:list") ){
				ifstream ppo_in(PRODART::ENV::get_option_value<string>("make_phi_psi_hist:list").c_str(), ios::in);

				if (!ppo_in.is_open()){
					cerr << "ERROR: can't open intput file: " << PRODART::ENV::get_option_value<string>("make_phi_psi_hist:list")
					                                        						   << endl;
					return -1;
				}
				else {
					cout << "making histogram from list" << endl;
					MISC::output_rama_histogram_from_raw_phi_psi_filelist(ppo_in,
							output_pdb_file,
							PRODART::ENV::get_option_value<bool>("make_phi_psi_hist:trans"));
				}

			}
			else {
				cerr << "ERROR: not enough options set"
						<< endl;
			}
		}
		else if (PRODART::ENV::get_option_value<bool>("make_dummy_pdb") == true
				&& output_pdb_is_set ){
			/*
			PRODART::ENV::register_option("make_dummy_pdb:seq",string("A"),"1 letter seq");
			PRODART::ENV::register_option("make_dummy_pdb:start_resnum",int(1),"start resnum");
			PRODART::ENV::register_option("make_dummy_pdb:dummy_chain",char('A'),"start resnum");
			*/
			pose_shared_ptr pose_ = PRODART::POSE::new_pose();
			pose_->add_new_chain(PRODART::ENV::get_option_value<char>("make_dummy_pdb:dummy_chain"));

			const string seq = PRODART::ENV::get_option_value<string>("make_dummy_pdb:seq");
			PRODART::POSE::residue_type_vector rtVec = POSE::one_letter_seq_to_residue_type_vector(seq);

			for (unsigned int i = 0; i < rtVec.size(); i++){
				residue_shared_ptr res = pose_->append_residue(rtVec[i], get_option_value<char>("make_dummy_pdb:dummy_chain"));
				pose_->index();
				const double rd = ENV::get_option_value<double>("make_dummy_pdb:rand");
				vector3d rand(ENV::get_random_num_gen()->rand(rd),ENV::get_random_num_gen()->rand(rd),ENV::get_random_num_gen()->rand(rd));
				pose_->add_new_atom(rand, atom_type("N"), res->get_internal_residue_index());
				rand = vector3d(ENV::get_random_num_gen()->rand(rd),ENV::get_random_num_gen()->rand(rd),ENV::get_random_num_gen()->rand(rd));
				pose_->add_new_atom(rand, atom_type("CA"), res->get_internal_residue_index());
				rand = vector3d(ENV::get_random_num_gen()->rand(rd),ENV::get_random_num_gen()->rand(rd),ENV::get_random_num_gen()->rand(rd));
				pose_->add_new_atom(rand, atom_type("C"), res->get_internal_residue_index());
				pose_->index();
			}

			pose_->renumber_residues(get_option_value<char>("make_dummy_pdb:dummy_chain"), get_option_value<int>("make_dummy_pdb:start_resnum"));
			pose_->index();
			pose_->outputPdb(output_pdb_file);

		}
		else if (ENV::get_option_value<bool>("make_centered_bb_frags") == true
				&& PRODART::ENV::is_set("pose:io:pdb:list:l")
				&& output_pdb_is_set){
			ifstream list_in(PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l").c_str(), ios::in);
			if (!list_in.is_open()){
				cerr << "ERROR: can't open input file list: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l")
			                                        						   << endl;
				return -1;
			}
			else {
				cout << "outputting fragments..." << endl;
				MISC::output_centred_backbone_segments(list_in,
						ENV::get_option_value<int>("make_centered_bb_frags:bb_frag_len"),
						output_pdb_file,
						ENV::get_option_value<bool>("make_centered_bb_frags:no_gly"),
						ENV::get_option_value<bool>("make_centered_bb_frags:no_pro"));
			}
		}
		else if (ENV::get_option_value<bool>("output_filtered_list") == true
				&& PRODART::ENV::is_set("pose:io:pdb:list:l")
				&& output_pdb_is_set){
			ifstream list_in(PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l").c_str(), ios::in);
			if (!list_in.is_open()){
				cerr << "ERROR: can't open input file list: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l")
			                                        						   << endl;
				return -1;
			}
			else {
				cout << "outputting filtered list..." << endl;
				MISC::output_quality_filtered_list(list_in, output_pdb_file);
			}
		}
		else if (ENV::get_option_value<bool>("process_pdbs_for_rosetta_match") == true
				&& PRODART::ENV::is_set("pose:io:pdb:list:l")
				&& PRODART::ENV::is_set("output_dir")){

			ifstream list_in(PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l").c_str(), ios::in);
			if (!list_in.is_open()){
				cerr << "ERROR: can't open input file list: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:list:l")
			                                        						   << endl;
				return -1;
			}
			else {
				cout << "processing files for Rosetta Match" << endl;
				MISC::process_list_for_rosetta_match(list_in, ENV::get_option_value<string>("output_dir"),
						ENV::get_option_value<int>("process_pdbs_for_rosetta_match:min_len"),
						ENV::get_option_value<int>("process_pdbs_for_rosetta_match:max_len"));
			}
		}
		else if (ENV::get_option_value<bool>("output_design_quality") == true){
			if (PRODART::ENV::is_set("output_design_quality:input")
				&& PRODART::ENV::is_set("output_design_quality:output") ){
				ifstream input(ENV::get_option_value<string>("output_design_quality:input").c_str());
				if (input.is_open()){
					ofstream output(ENV::get_option_value<string>("output_design_quality:output").c_str());
					if (output.is_open()){
						MISC::output_design_quality(input, output);
					}
					else {
						cout << "ERROR: couldn't open --output_design_quality:input file: " << ENV::get_option_value<string>("output_design_quality:input") << endl;
					}
				}
				else {
					cout << "ERROR: couldn't open --output_design_quality:output file: " << ENV::get_option_value<string>("output_design_quality:output") << endl;
				}
			}
			else {
				cerr << "ERROR: you need to set --output_design_quality:input and --output_design_quality:output with this option" << endl;
			}
		}
		else {
			cerr << "Nothing to do it's up to you" << endl;
		}
	}

	return 0;
}

