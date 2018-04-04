/*
 * pd2_ca2main.cpp
 *
 *  Created on: 24 Sep 2010
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
#include "backbonebuilder/backbone_builder.h"
#include "protocols/ca2main_protocols.h"
#include "protocols/misc_protocols.h"

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

using namespace PRODART::ENV;
using namespace PRODART::UTILS;
using namespace boost;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;
using namespace PRODART::PROTOCOLS;
using namespace std;


typedef std::vector< std::string >  split_vector_type;


using PRODART::UTILS::PI;
using PRODART::POSE::atom_type;
using PRODART::UTILS::vector3d;

using namespace PRODART::POSE;
using namespace PRODART;

int main( int argc, char *argv[] ) {



	//options_manager &opt_mgr = prodart_env::Instance()->get_options_manager();
	//PRODART::ENV::list_options_full_format(cout);
	PRODART::ENV::register_option("help",bool(false),"print options");
	PRODART::ENV::register_option("ca2main:simple",bool(false),"simple ca2main with no CA refinement but including BB minimisation");
	PRODART::ENV::register_option("ca2main:simple_secs_refine",bool(false),"simple ca2main including CA refinement and BB minimisation");
	PRODART::ENV::register_option("ca2main:simple_secs_indep_refine",bool(false),"simple ca2main including CA refinement and BB minimisation with independent sse refinement");
	PRODART::ENV::register_option("ca2main:alphabet_fixed_ca",
			bool(false),"alphabet ca2main with no CA refinement but including BB minimisation (excl CA)", true);

	PRODART::ENV::register_option("ca2main:new",bool(false),"new alphabet ca2main with no CA refinement but including BB minimisation");
	PRODART::ENV::register_option("ca2main:new_fixed_ca",bool(false),"new alphabet ca2main with no CA refinement but including BB minimisation");


	PRODART::ENV::register_option("ca2main:secs",string(""),"sec struct file");
	PRODART::ENV::register_option("ca2main:bb_min_steps",(unsigned long)(500),"backbone min steps");
	PRODART::ENV::register_option("ca2main:ca_mc_steps",(unsigned long)(100000),"ca monte carlo refinements steps");


	PRODART::ENV::register_option("ca2main:alphabet",string(""),"alphabet input");

	PRODART::ENV::register_option("ca2main:chain",char(' '),"chain ID for residue range");
	PRODART::ENV::register_option("ca2main:start_resnum",string(""),"start resnum for residue range");
	PRODART::ENV::register_option("ca2main:end_resnum",string(""),"end resnum for residue range");

	PRODART::ENV::register_option("ca2main:force_old",bool(false),"FOR DEBUGGING only: force old ca2main method");



	ENV::unhide_option("restraints::rstfile");

	if (!prodart_env::Instance()->init(argc, argv)){
		cerr << "initialisation errors" << endl;
		return -1;
	}

	cout << PRODART::ENV::get_full_cmd_line() << endl;

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

	PRODART::POSE::three_state_sec_struct_vector secs;
	if ( PRODART::ENV::is_set("ca2main:secs")) {
		secs = MISC::read_3state_secs(PRODART::ENV::get_option_value<string>("ca2main:secs"));
		cout << secs << endl;
	}
	PRODART::POSE::pose_shared_ptr in_pdb = PRODART::POSE::new_pose();
	ifstream protein_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:i").c_str(), ios::in);
	if (input_pdb_is_set){
		if (protein_file.is_open()){
			in_pdb->loadPdb(protein_file);
			in_pdb->set_label(PRODART::ENV::get_option_value<string>("pose:io:pdb:i"));
		}
		else {
			cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:i")
			                                           << endl;
			return -1;
		}
		protein_file.close();
	}
	else{
		cerr << "ERROR: no pdb input file set\n"
			 << endl;
		return -1;
	}


	ofstream output_pdb_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:o").c_str(), ios::out);
	if (output_pdb_is_set){
		if (!output_pdb_file.is_open()){
			cerr << "ERROR: can't open output file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:o")
			                                           << endl;
			return -1;
		}
	}
	else {
		cerr << "ERROR: no pdb input file set\n"
			 << endl;
		return -1;
	}

	if (PRODART::ENV::get_option_value<bool>("ca2main:simple") == true){
		validate_ca_pose(in_pdb);
		if (CA2MAIN::simple_ca2main_minimise(in_pdb, ENV::get_option_value<unsigned long>("ca2main:bb_min_steps"))){
			cout << "Validating output PDB:" << endl;
			const bool valid = validate_bb_pose(in_pdb);
			if (!valid){
				cout << "WARNING: output PDB has bad bonds:" << endl;
			}
			else {
				cout << "output PDB seems OK" << endl;
			}
			in_pdb->outputPdb(output_pdb_file);
		}
		else {
			cerr << "ERROR: couldn't add mainchain atoms\n"
				 << endl;
			return -1;
		}

		return 0;

	}
	if (PRODART::ENV::get_option_value<bool>("ca2main:simple_secs_refine") == true){
		if (!PRODART::ENV::is_set("ca2main:secs")){
			cerr << "ERROR: --ca2main:secs option not set" << endl;
			return -1;
		}
		validate_ca_pose(in_pdb);
		//cout << "not implemented yet" << endl;
		if (CA2MAIN::simple_ca2main_secs_refine_minimise(in_pdb, secs,
				ENV::get_option_value<unsigned long>("ca2main:ca_mc_steps"),
				ENV::get_option_value<unsigned long>("ca2main:bb_min_steps"))){
			cout << "Validating output PDB:" << endl;
			const bool valid = validate_bb_pose(in_pdb);
			if (!valid){
				cout << "WARNING: output PDB has bad bonds:" << endl;
			}
			else {
				cout << "output PDB seems OK" << endl;
			}
			in_pdb->outputPdb(output_pdb_file);
		}
		else {
			cerr << "ERROR: couldn't add mainchain atoms\n"
				 << endl;
			return -1;
		}

		return 0;
	}
	if (PRODART::ENV::get_option_value<bool>("ca2main:simple_secs_indep_refine") == true){
		if (!PRODART::ENV::is_set("ca2main:secs")){
			cerr << "ERROR: --ca2main:secs option not set" << endl;
			return -1;
		}
		validate_ca_pose(in_pdb);
		//cout << "not implemented yet" << endl;
		if (CA2MAIN::simple_ca2main_secs_indep_refine_minimise(in_pdb, secs)){
			cout << "Validating output PDB:" << endl;
			const bool valid = validate_bb_pose(in_pdb);
			if (!valid){
				cout << "WARNING: output PDB has bad bonds:" << endl;
			}
			else {
				cout << "output PDB seems OK" << endl;
			}
			in_pdb->outputPdb(output_pdb_file);
		}
		else {
			cerr << "ERROR: couldn't add mainchain atoms\n"
				 << endl;
			return -1;
		}

		return 0;
	}
	if (PRODART::ENV::get_option_value<bool>("ca2main:alphabet_fixed_ca") == true){
		if( PRODART::ENV::is_set("ca2main:alphabet")){
			ifstream fbb_aleph(ENV::get_option_value<string>("ca2main:alphabet").c_str(), ios::in);
			BB_BUILDER::alphabet_bb_builder_shared_ptr bb_aleph = BB_BUILDER::new_alphabet_bb_builder();
			if (fbb_aleph.is_open()){
				bb_aleph->load_alphabet_pose(fbb_aleph);
			}
			else {
				cerr << "ERROR: could not open alphabet file: " << ENV::get_option_value<string>("ca2main:alphabet") << endl;
			}
			validate_ca_pose(in_pdb);
			if (CA2MAIN::alphabet_ca2main_fixed_ca_minimise(in_pdb, bb_aleph, ENV::get_option_value<unsigned long>("ca2main:bb_min_steps"))){
				cout << "Validating output PDB:" << endl;
				const bool valid = validate_bb_pose(in_pdb);
				if (!valid){
					cout << "WARNING: output PDB has bad bonds:" << endl;
				}
				else {
					cout << "output PDB seems OK" << endl;
				}
				in_pdb->outputPdb(output_pdb_file);
			}
			else {
				cerr << "ERROR: couldn't add mainchain atoms\n"
						<< endl;
				return -1;
			}

			return 0;
		}
		else {
			cerr << "ERROR: option --ca2main:alphabet required" << endl;
		}

	}
	if (PRODART::ENV::get_option_value<bool>("ca2main:new") == true){
		validate_ca_pose(in_pdb);
		if (CA2MAIN::new_ca2main_minimise(in_pdb, ENV::get_option_value<unsigned long>("ca2main:bb_min_steps"))){
			cout << "Validating output PDB:" << endl;
			const bool valid = validate_bb_pose(in_pdb);
			if (!valid){
				cout << "WARNING: output PDB has bad bonds:" << endl;
			}
			else {
				cout << "output PDB seems OK" << endl;
			}
			in_pdb->outputPdb(output_pdb_file);
		}
		else {
			cerr << "ERROR: couldn't add mainchain atoms\n"
				 << endl;
			return -1;
		}

		return 0;

	}
	if (PRODART::ENV::get_option_value<bool>("ca2main:new_fixed_ca") == true){
		validate_ca_pose(in_pdb);
		if (!ENV::is_set("ca2main:start_resnum")
			&& !ENV::is_set("ca2main:end_resnum")){
			if (CA2MAIN::new_ca2main_minimise_fixed_ca(in_pdb, ENV::get_option_value<unsigned long>("ca2main:bb_min_steps"))){
				cout << "Validating output PDB:" << endl;
				const bool valid = validate_bb_pose(in_pdb);
				if (!valid){
					cout << "WARNING: output PDB has bad bonds:" << endl;
				}
				else {
					cout << "output PDB seems OK" << endl;
				}
				in_pdb->outputPdb(output_pdb_file);
			}
			else {
				cerr << "ERROR: couldn't add mainchain atoms\n"
						<< endl;
				return -1;
			}
		}
		else if (ENV::is_set("ca2main:start_resnum")
			&& ENV::is_set("ca2main:end_resnum")){
			cout << "running ca2main:new_fixed_ca with residue range\n"
					<< endl;
			residue_shared_ptr start_res = in_pdb->get_residue(ENV::get_option_value<string>("ca2main:start_resnum"),
					ENV::get_option_value<char>("ca2main:chain"));
			residue_shared_ptr end_res = in_pdb->get_residue(ENV::get_option_value<string>("ca2main:end_resnum"),
					ENV::get_option_value<char>("ca2main:chain"));

			bool_vector loop_mask = PROTOCOLS::CA2MAIN::make_residue_loop_mask(in_pdb,
					start_res->get_internal_residue_index(),
					end_res->get_internal_residue_index());
			//PRINT_EXPR(loop_mask);
			for (unsigned int ii = 0 ; ii < loop_mask.size(); ii++){
				cout << loop_mask[ii];
			}
			cout << endl;
			for (unsigned int ii = 0 ; ii < loop_mask.size(); ii++){
				cout << (ii % 10);
			}
			cout << endl;
			if (!ENV::get_option_value<bool>("ca2main:force_old") ? CA2MAIN::new_ca2main_minimise_fixed_ca(in_pdb,
					loop_mask,
					ENV::get_option_value<unsigned long>("ca2main:bb_min_steps"))
					: CA2MAIN::simple_ca2main_minimise_fixed_ca(in_pdb,
							loop_mask,
							ENV::get_option_value<unsigned long>("ca2main:bb_min_steps"))){
				/*
				cout << "Validating output PDB:" << endl;
				const bool valid = validate_bb_pose(in_pdb);
				if (!valid){
					cout << "WARNING: output PDB has bad bonds:" << endl;
				}
				else {
					cout << "output PDB seems OK" << endl;
				}
				*/
				in_pdb->outputPdb(output_pdb_file);
			}
			else {
				cerr << "ERROR: couldn't add mainchain atoms\n"
						<< endl;
				return -1;
			}

		}
		else {
			cerr << "ERROR: both start_resnum and end_resnum need to be set\n"
				 << endl;
			return -1;
		}

		return 0;

	}
	else {
		cerr << "Nothing to do it's up to you" << endl;
	}

	return 0;

}






