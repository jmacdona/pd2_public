/*
 * pd2_mpi_monte_carlo.cpp
 *
 *  Created on: 29 Sep 2013
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
	PRODART::ENV::register_option("mpi_mc:ca_mc",bool(false),"mc");
	PRODART::ENV::register_option("mpi_mc:steps",(unsigned long)(10000),"minimisation steps");
	PRODART::ENV::register_option("mpi_mc:ca_pot",string("ca_default"),"ca_potential name", true);
	PRODART::ENV::register_option("mpi_mc:beta",double(1.0),"beta");
	PRODART::ENV::register_option("mpi_mc:replica_specs",string(""),"beta");
	PRODART::ENV::register_option("mpi_mc:exchange_freq",(unsigned long)(4000),"exchange frequency");
	PRODART::ENV::register_option("mpi_mc:start_extended",bool(false),"start from extended conformation");



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

	PRODART::POSE::pose_shared_ptr ref_pdb = in_pdb->clone();



	if (ENV::get_option_value<bool>("mpi_mc:start_extended") == true){
		cout << "Extending chain" << endl;
		MISC::extend_ca_chain(in_pdb);
		cout << "chain extended" << endl;

		{
			string output_pdb_filename = PRODART::ENV::get_option_value<string>("output_root");
			output_pdb_filename.append(".replica");
			output_pdb_filename.append(boost::lexical_cast<string>(world.rank()));
			output_pdb_filename.append(".start.pdb");

			ofstream output_pdb_file(output_pdb_filename.c_str(), ios::out);
			if (!output_pdb_file.is_open()){
				cerr << "ERROR: can't open output file: " << output_pdb_filename
						<< endl;
				return -1;
			}
			cout << "outputting to: " << output_pdb_filename << endl;
			in_pdb->outputPdb(output_pdb_file, true);
		}
	}



	if (PRODART::ENV::get_option_value<bool>("mpi_mc:ca_mc") == true){

		potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot = pot_factory->make_preset_potentials_container(PRODART::ENV::get_option_value<string>("mpi_mc:ca_pot"));

		PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(in_pdb);

		PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = MOVERS::move_set_factory::Instance()->make_preset_move_set(meta_data,
				"ca_standard",
				PRODART::ENV::get_random_num_gen());

		PRODART::POSE::SIM::MPI::mpi_replica_exchange_mc_protocol_verbose_shared_ptr protocol = PRODART::POSE::SIM::MPI::new_mpi_replica_exchange_mc_protocol_verbose( world,
				ENV::get_option_value<unsigned long>("mpi_mc:steps"),
				ENV::get_random_num_gen());

		protocol->set_exchange_frequency(ENV::get_option_value<unsigned long>("mpi_mc:exchange_freq"));

		PRODART::POSE::SIM::monte_carlo_shared_ptr mc_sim = PRODART::POSE::SIM::new_monte_carlo(pot,
				overall_move_set,
				protocol);

		if (protocol->load_replica_specs(ENV::get_option_value<string>("mpi_mc:replica_specs"))){

			cout << "running ca anneal" << endl;

			mc_sim->make_move(meta_data);

			const int replica_number = mc_sim->get_state().replica_number;

			string output_pdb_filename = PRODART::ENV::get_option_value<string>("output_root");
			output_pdb_filename.append(".replica");
			output_pdb_filename.append(boost::lexical_cast<string>(replica_number));
			output_pdb_filename.append(".pdb");

			ofstream output_pdb_file(output_pdb_filename.c_str(), ios::out);
			if (!output_pdb_file.is_open()){
				cerr << "ERROR: can't open output file: " << output_pdb_filename
						<< endl;
				return -1;
			}

			in_pdb->outputPdb(output_pdb_file, true);

			mc_sim->get_state().print_summary(cout);
		}
		else {
			cerr << "ERROR: some problem with replica specs file or mpirun arguments: " << ENV::get_option_value<string>("mpi_mc:replica_specs")
					<< endl;
		}

	}
	else {
		cerr << "Nothing to do it's up to you" << endl;
	}

	return 0;

}
