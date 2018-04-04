/*
 * pd2_minimise.cpp
 *
 *  Created on: Nov 16, 2010
 *      Author: jmacdon
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
#include "protocols/loop_model_protocols.h"


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
	PRODART::ENV::register_option("minimise:minimise",bool(false),"simple minimisation");
	PRODART::ENV::register_option("minimise:steps",(unsigned long)(500),"minimisation steps");
	PRODART::ENV::register_option("minimise:bb_pot",string("bb_min_default"),"bb_potential name");


	PRODART::ENV::register_option("minimise:ca_minimise",bool(false),"simple minimisation");

	PRODART::ENV::unhide_option("loop_model:loop_mask");


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



	if (PRODART::ENV::get_option_value<bool>("minimise:minimise") == true && ENV::is_set("loop_model:loop_mask") == false){
		cout << "pd2_minimise: minimising whole backbone..." << endl;
		POSE_UTILS::quick_add_HN(in_pdb, false);

		if (CA2MAIN::bb_minimise(in_pdb, PRODART::ENV::get_option_value<string>("minimise:bb_pot"), PRODART::ENV::get_option_value<unsigned long>("minimise:steps"))){	//bb_minimise_geom(in_pdb, PRODART::ENV::get_option_value<unsigned long>("minimise:steps"))){
			in_pdb->outputPdb(output_pdb_file);
			if (!validate_bb_pose(in_pdb)){
				cerr << "pd2_minimise: WARNING: could not validate bb geometry" << endl;
			}
			else {
				cout << "pd2_minimise: bb geometry seems OK" << endl;
			}
		}
		else {
			cerr << "ERROR: couldn't add mainchain atoms\n"
				 << endl;
			return -1;
		}

		return 0;

	}
	else if (PRODART::ENV::get_option_value<bool>("minimise:minimise") == true && ENV::is_set("loop_model:loop_mask") == true){
		cout << "pd2_minimise: minimising loop backbone only..." << endl;
		const PROTOCOLS::LOOP::bool_vector loop_mask = PROTOCOLS::MISC::make_bool_vec(PRODART::ENV::get_option_value<string>("loop_model:loop_mask"));
		POSE_UTILS::quick_add_HN(in_pdb, false);
		if (CA2MAIN::bb_minimise(in_pdb, loop_mask, PRODART::ENV::get_option_value<string>("minimise:bb_pot"), PRODART::ENV::get_option_value<unsigned long>("minimise:steps"))){	//bb_minimise_geom(in_pdb, PRODART::ENV::get_option_value<unsigned long>("minimise:steps"))){
			in_pdb->outputPdb(output_pdb_file);
			if (!validate_bb_pose(in_pdb)){
				cerr << "pd2_minimise: WARNING: could not validate bb geometry" << endl;
			}
			else {
				cout << "pd2_minimise: bb geometry seems OK" << endl;
			}
		}
		else {
			cerr << "ERROR: couldn't add mainchain atoms\n"
				 << endl;
			return -1;
		}

		return 0;

	}
	else if (PRODART::ENV::get_option_value<bool>("minimise:ca_minimise") == true && ENV::is_set("loop_model:loop_mask") == false){
		cout << "pd2_minimise: minimising ca bonds and bumps..." << endl;

		if (LOOP::ca_minimise_bonds_bumps(in_pdb)){	//bb_minimise_geom(in_pdb, PRODART::ENV::get_option_value<unsigned long>("minimise:steps"))){
			ca_pose_meta::inactivate_pseudo_NO(in_pdb);
			in_pdb->outputPdb(output_pdb_file);
			if (!validate_ca_pose(in_pdb)){
				cerr << "pd2_minimise: WARNING: could not validate ca geometry" << endl;
			}
			else {
				cout << "pd2_minimise: bb geometry seems OK" << endl;
			}
		}
		else {
			cerr << "ERROR: couldn't add mainchain atoms\n"
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
