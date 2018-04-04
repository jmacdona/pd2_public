/*
 * pd2_mpi_conpred.cpp
 *
 *  Created on: 9 Oct 2013
 *      Author: jmacdona
 */



#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>


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
#include "movers/move_set_factory.h"
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
#include "simulation/mpi/mpi_replica_exchange_mc_protocol.h"
#include "simulation/mpi/mpi_replica_exchange_mc_protocol_verbose.h"
#include "backbonebuilder/backbone_builder.h"
#include "protocols/ca2main_protocols.h"
#include "protocols/misc_protocols.h"
#include "protocols/mpi/structure_prediction.h"

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
	mpi::environment env(argc, argv);
	mpi::communicator world;


	//options_manager &opt_mgr = prodart_env::Instance()->get_options_manager();
	//PRODART::ENV::list_options_full_format(cout);
	PRODART::ENV::register_option("help",bool(false),"print options");
	PRODART::ENV::register_option("conpred:remc_steps",(unsigned long)(10000),"replica mc steps steps");
	PRODART::ENV::register_option("conpred:replica_specs",string(""),"beta");
	//PRODART::ENV::register_option("conpred:exchange_freq",(unsigned long)(4000),"exchange frequency");
	PRODART::ENV::register_option("conpred:secs",string(""),"sec struct file");
	PRODART::ENV::register_option("conpred:seq",string(),"protein sequence");
	PRODART::ENV::register_option("conpred:psicov_in",string(),"psicov contact predictions");
	PRODART::ENV::hide_option("pose:io:pdb:i");
	PRODART::ENV::hide_option("pose:io:pdb:o");
	//PRODART::ENV::register_option("mpi_mc:replica_specs",string(""),"beta");

	if (!prodart_env::Instance()->init(argc, argv)){
		cerr << "initialisation errors" << endl;
		return -1;
	}



	if (world.rank() == 0)	cout << PRODART::ENV::get_full_cmd_line() << endl;
		if (PRODART::ENV::get_option_value<bool>("help") == true){
			if (world.rank() == 0) PRODART::ENV::list_options_full_format(cout);
			return 0;
		}

	if (PRODART::ENV::is_set("")){
		if (world.rank() == 0) cerr << "ERROR: command line free options are not supported " << world.rank() << endl;
		return -1;
	}



	POSE::pose_shared_ptr pose_ = new_pose();
	PRODART::POSE::residue_type_vector seq;
	if (ENV::is_set("conpred:seq") && !ENV::is_set("pose:io:pdb:i")){
		seq = POSE::one_letter_seq_to_residue_type_vector(ENV::get_option_value<string>("conpred:seq"));
	}
	else if (ENV::is_set("pose:io:pdb:i") && !ENV::is_set("conpred:seq")) {
		ifstream protein_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:i").c_str(), ios::in);
		if (protein_file.is_open()){
			pose_->loadPdb(protein_file);
			pose_->set_label(PRODART::ENV::get_option_value<string>("pose:io:pdb:i"));
			if (pose_->get_chain_count() ==1){
				seq = pose_->get_chain(0)->get_sequence();
			}
			else {
				cerr << "ERROR: can only deal with single chain structures" << endl;
				return -1;
			}
		}
		else {
			cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:i")
			                                           << endl;
			return -1;
		}
	}
	else {
		// ERROR
		cerr << "ERROR: need sequence from either --conpred:secs or --pose:io:pdb:i options but not both" << endl;
		return -1;
	}

	PRODART::POSE::three_state_sec_struct_vector secs;
	if ( PRODART::ENV::is_set("conpred:secs")) {
		secs = MISC::read_3state_secs(PRODART::ENV::get_option_value<string>("conpred:secs"));
		cout << secs << endl;
	}
	else {
		cerr << "WARNING: secondary structure prediction not set with --secs option" << endl;
		secs = PRODART::POSE::three_state_sec_struct_vector(seq.size(), ss3_OTHER);
	}

	if (!ENV::is_set("conpred:psicov_in")){
		cerr << "ERROR: need PSICOV predicted contacts set --conpred:psicov_in" << endl;
		return -1;
	}


	pose_shared_ptr start_struct;
	if (ENV::is_set("pose:io:pdb:i")){
		start_struct = pose_;
	}



	if (PROTOCOLS::MPI::STRUCTURE_PREDICTION::mpi_replica_exchange_folding(seq,
			secs,
			ENV::get_option_value<string>("conpred:psicov_in"),
			ENV::get_option_value<unsigned long>("conpred:remc_steps"),
			world,
			ENV::get_option_value<string>("conpred:replica_specs"),
			start_struct)){

		cout << "finished... " << endl;

	}
	else {
		cerr << "ERROR: " << endl;
	}




	return 0;

}




