/*
 * pd2_loop_model.cpp
 *
 *  Created on: Oct 4, 2010
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


void print_phi_psi(pose_shared_ptr pose_, int start, int end, ostream& output );

int main( int argc, char *argv[] ) {



	//options_manager &opt_mgr = prodart_env::Instance()->get_options_manager();
	//PRODART::ENV::list_options_full_format(cout);
	PRODART::ENV::register_option("help",bool(false),"print options");
	PRODART::ENV::register_option("loop_model:ca_insert",bool(false),"insert loop");
	PRODART::ENV::register_option("loop_model:fiser_test_protocol",bool(false),"fiser test protocol -DO NOT USE - USE: fiser_test_protocol_corrected");
	PRODART::ENV::register_option("loop_model:fiser_test_protocol_corrected",bool(false),"fiser test protocol");
	PRODART::ENV::register_option("loop_model:fiser_test_protocol_control",bool(false),"fiser test protocol control - DO NOT USE - USE: fiser_test_protocol_control_corrected");
	PRODART::ENV::register_option("loop_model:fiser_test_protocol_control_corrected",bool(false),"fiser test protocol control");
	PRODART::ENV::register_option("loop_model:fiser_print_phi_psi_omega",bool(false),"print phi psi omega angles");
	PRODART::ENV::register_option("loop_model:fiser_print_phi_psi_omega_corrected",bool(false),"print phi psi omega angles");
	PRODART::ENV::register_option("loop_model:filter_loops", bool(true),"filter loops during loop sampling");
	PRODART::ENV::register_option("loop_model:fiser_bb_loop_rmsd",bool(false),"get loop global RMSD");
	PRODART::ENV::register_option("loop_model:fiser_bb_loop_rmsd_corrected",bool(false),"get loop global RMSD");

	PRODART::ENV::register_option("loop_model:insert_seq",string(),"insert sequence");
	PRODART::ENV::register_option("loop_model:insert_overwrite",int(0),"number to overwrite");
	PRODART::ENV::register_option("loop_model:insert_resnum",string(""),"insert position residue number (inserts before this residue)");
	PRODART::ENV::register_option("loop_model:insert_chain",char(' '),"insert position chain ID (inserts before this residue)");
	PRODART::ENV::register_option("loop_model:num_structs",(unsigned long)(1),"number of structures to output");
	PRODART::ENV::register_option("loop_model:renumber",bool(false),"renumber?");
	PRODART::ENV::register_option("loop_model:seq_dep",bool(false),"sequence dependence");
	ENV::register_option("loop_model:bb_min_steps",(unsigned long)(100),"number of bb min steps");



	if (!prodart_env::Instance()->init(argc, argv)){
		cerr << "initialisation errors" << endl;
		return -1;
	}

	cout << PRODART::ENV::get_full_cmd_line() << endl;

	if (PRODART::ENV::get_option_value<bool>("help") == true){
		PRODART::ENV::list_options_full_format(cout);
		return 0;
	}

	//cout << PRODART::ENV::get_full_cmd_line() << endl;


	if (PRODART::ENV::is_set("")){
		cerr << "ERROR: command line free options are not supported" << endl;
		return -1;
	}

	PRODART::ENV::list_options_full_format(cout);

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
		//protein_file.close();
	}
	else{
		cerr << "ERROR: no pdb input file set\n"
			 << endl;
		return -1;
	}

	ofstream output_pdb_file;
	if (!(PRODART::ENV::get_option_value<bool>("loop_model:fiser_bb_loop_rmsd") || PRODART::ENV::get_option_value<bool>("loop_model:fiser_bb_loop_rmsd_corrected"))){
		output_pdb_file.open(PRODART::ENV::get_option_value<string>("pose:io:pdb:o").c_str(), ios::out);
		if (output_pdb_is_set){
			if (!output_pdb_file.is_open()){
				cerr << "ERROR: can't open output file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:o")
			                                        		   << endl;
				return -1;
			}
		}
		else {
			cerr << "ERROR: no pdb output file set\n"
					<< endl;
			return -1;
		}
	}



	if (PRODART::ENV::get_option_value<bool>("loop_model:ca_insert") == true){

		const string seq = PRODART::ENV::get_option_value<string>("loop_model:insert_seq");

		PRODART::POSE::residue_type_vector rtVec = POSE::one_letter_seq_to_residue_type_vector(seq);

		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const string ca_pot_str = PRODART::ENV::get_option_value<bool>("loop_model:seq_dep") ? string("ca_seq_dep_default") : string("ca_default");


			PROTOCOLS::LOOP::remodel_loops_multi_filtered(in_pdb,
					res->get_internal_residue_index(),
					rtVec,
					overwrite, //4,
					50,		//50
					100,	//300
					250,	//400
					0.2,
					1.2,
					PRODART::ENV::get_option_value<unsigned long>("loop_model:num_structs"),
					output_pdb_file,
					PRODART::ENV::get_option_value<bool>("loop_model:renumber"),
					ENV::get_option_value<unsigned long>("loop_model:bb_min_steps"),
					ca_pot_str);

			/*
			bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(in_pdb);
			cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
			cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
			*/

			//in_pdb->outputPdb(output_pdb_file);
			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}


	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_test_protocol") == true){
		if (PRODART::ENV::is_set("loop_model:insert_seq")){
			cout << "option --loop_model:insert_seq is ignored" << endl;
		}



		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();




			PRODART::POSE::residue_type_vector rtVec;// = POSE::one_letter_seq_to_residue_type_vector(seq);

			for (int i = resnum; i < resnum + overwrite; i++){
				rtVec.push_back(in_pdb->get_residue(i)->get_type());
			}

			string loop_seq = POSE::residue_type_map::Instance()->to_1_letter_code(rtVec);

			cout << "loop_sequence:\t" << loop_seq << endl ;

			const string ca_pot_str = PRODART::ENV::get_option_value<bool>("loop_model:seq_dep") ? string("ca_seq_dep_default") : string("ca_default");

			PROTOCOLS::LOOP::remodel_loops_multi_filtered(in_pdb,
					res->get_internal_residue_index(),
					rtVec,
					overwrite, //4,
					50,		//50
					100,	//300
					250,	//400
					0.2,
					1.2,
					PRODART::ENV::get_option_value<unsigned long>("loop_model:num_structs"),
					output_pdb_file,
					PRODART::ENV::get_option_value<bool>("loop_model:renumber"),
					ENV::get_option_value<unsigned long>("loop_model:bb_min_steps"),
					ca_pot_str,
					"bb_min_default",
					!PRODART::ENV::get_option_value<bool>("loop_model:filter_loops"));

			/*
			bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(in_pdb);
			cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
			cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
			*/

			//in_pdb->outputPdb(output_pdb_file);
			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}


	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_test_protocol_corrected") == true){
		if (PRODART::ENV::is_set("loop_model:insert_seq")){
			cout << "option --loop_model:insert_seq is ignored" << endl;
		}



		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		if (overwrite < 3){
			cerr << "ERROR: loop_model:fiser_test_protocol_corrected: minimum loop size is 3" << endl;
			return -1;
		}

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();




			PRODART::POSE::residue_type_vector rtVec;// = POSE::one_letter_seq_to_residue_type_vector(seq);

			for (int i = resnum+1; i < resnum + overwrite-1; i++){
				rtVec.push_back(in_pdb->get_residue(i)->get_type());
			}

			string loop_seq = POSE::residue_type_map::Instance()->to_1_letter_code(rtVec);

			cout << "loop_sequence:\t" << loop_seq << endl ;

			const string ca_pot_str = PRODART::ENV::get_option_value<bool>("loop_model:seq_dep") ? string("ca_seq_dep_default") : string("ca_default");

			PROTOCOLS::LOOP::remodel_loops_multi_filtered(in_pdb,
					res->get_internal_residue_index() + 1,
					rtVec,
					overwrite - 2, //4,
					50,		//50
					100,	//300
					250,	//400
					0.2,
					1.2,
					PRODART::ENV::get_option_value<unsigned long>("loop_model:num_structs"),
					output_pdb_file,
					PRODART::ENV::get_option_value<bool>("loop_model:renumber"),
					ENV::get_option_value<unsigned long>("loop_model:bb_min_steps"),
					ca_pot_str,
					"bb_min_default",
					!PRODART::ENV::get_option_value<bool>("loop_model:filter_loops"));

			/*
			bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(in_pdb);
			cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
			cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
			*/

			//in_pdb->outputPdb(output_pdb_file);
			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}


	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_test_protocol_control") == true){
		if (PRODART::ENV::is_set("loop_model:insert_seq")){
			cout << "option --loop_model:insert_seq is ignored" << endl;
		}



		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();




			PRODART::POSE::residue_type_vector rtVec;// = POSE::one_letter_seq_to_residue_type_vector(seq);

			for (int i = resnum; i < resnum + overwrite; i++){
				rtVec.push_back(in_pdb->get_residue(i)->get_type());
			}

			string loop_seq = POSE::residue_type_map::Instance()->to_1_letter_code(rtVec);

			cout << "loop_sequence:\t" << loop_seq << endl ;

			const string ca_pot_str = "ca_bond_bump";//PRODART::ENV::get_option_value<bool>("loop_model:seq_dep") ? string("ca_seq_dep_default") : string("ca_default");

			PROTOCOLS::LOOP::remodel_loops_multi_filtered(in_pdb,
					res->get_internal_residue_index(),
					rtVec,
					overwrite, //4,
					50,		//50
					100,	//300
					250,	//400
					0.2,
					1.2,
					PRODART::ENV::get_option_value<unsigned long>("loop_model:num_structs"),
					output_pdb_file,
					PRODART::ENV::get_option_value<bool>("loop_model:renumber"),
					ENV::get_option_value<unsigned long>("loop_model:bb_min_steps"),
					ca_pot_str,
					"bb_min_default",
					!PRODART::ENV::get_option_value<bool>("loop_model:filter_loops"));

			/*
			bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(in_pdb);
			cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
			cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
			*/

			//in_pdb->outputPdb(output_pdb_file);
			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}


	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_test_protocol_control_corrected") == true){
		if (PRODART::ENV::is_set("loop_model:insert_seq")){
			cout << "option --loop_model:insert_seq is ignored" << endl;
		}



		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		if (overwrite < 3){
			cerr << "ERROR: loop_model:fiser_test_protocol_corrected: minimum loop size is 3" << endl;
			return -1;
		}

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();




			PRODART::POSE::residue_type_vector rtVec;// = POSE::one_letter_seq_to_residue_type_vector(seq);

			for (int i = resnum+1; i < resnum + overwrite-1; i++){
				rtVec.push_back(in_pdb->get_residue(i)->get_type());
			}

			string loop_seq = POSE::residue_type_map::Instance()->to_1_letter_code(rtVec);

			cout << "loop_sequence:\t" << loop_seq << endl ;

			const string ca_pot_str = "ca_bond_bump";//PRODART::ENV::get_option_value<bool>("loop_model:seq_dep") ? string("ca_seq_dep_default") : string("ca_default");

			PROTOCOLS::LOOP::remodel_loops_multi_filtered(in_pdb,
					res->get_internal_residue_index(),
					rtVec,
					overwrite-2, //4,
					50,		//50
					100,	//300
					250,	//400
					0.2,
					1.2,
					PRODART::ENV::get_option_value<unsigned long>("loop_model:num_structs"),
					output_pdb_file,
					PRODART::ENV::get_option_value<bool>("loop_model:renumber"),
					ENV::get_option_value<unsigned long>("loop_model:bb_min_steps"),
					ca_pot_str,
					"bb_min_default",
					!PRODART::ENV::get_option_value<bool>("loop_model:filter_loops"));

			/*
			bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(in_pdb);
			cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
			cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
			*/

			//in_pdb->outputPdb(output_pdb_file);
			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}


	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_print_phi_psi_omega") == true){


		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();



			print_phi_psi(in_pdb, resnum-1, resnum + overwrite, output_pdb_file);
			while (!protein_file.eof()) {
				pose_shared_ptr pose_ = PRODART::POSE::new_pose();;
				pose_->loadPdb(protein_file);
				if (pose_->get_residue_count() > resnum + overwrite){
					print_phi_psi(pose_, resnum-1, resnum + overwrite, output_pdb_file);
				}
			}


			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}
	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_print_phi_psi_omega_corrected") == true){


		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();



			print_phi_psi(in_pdb, resnum, resnum + overwrite - 1, output_pdb_file);
			while (!protein_file.eof()) {
				pose_shared_ptr pose_ = PRODART::POSE::new_pose();;
				pose_->loadPdb(protein_file);
				if (pose_->get_residue_count() > resnum + overwrite){
					print_phi_psi(pose_, resnum, resnum + overwrite - 1, output_pdb_file);
				}
			}


			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}
	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_bb_loop_rmsd") == true){

		if (!ENV::is_set("pose:io:pdb:r")){
			cerr << "ERROR: reference pdb structure not set\n";
			return -1;
		}

		PRODART::POSE::pose_shared_ptr ref_pdb = PRODART::POSE::new_pose();
		ifstream ref_protein_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:r").c_str(), ios::in);

		if (ref_protein_file.is_open()){
			ref_pdb->loadPdb(ref_protein_file);
			ref_pdb->set_label(PRODART::ENV::get_option_value<string>("pose:io:pdb:r"));
		}
		else {
			cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:r")
				                                        		   << endl;
			return -1;
		}
		ref_protein_file.close();



		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();

			const double rmsd = PROTOCOLS::LOOP::get_fiser_bb_loop_rmsd(ref_pdb, in_pdb, resnum, overwrite);

			cout << "LOOP_RMSD:\t" << rmsd << endl;

			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}
	}
	else if (PRODART::ENV::get_option_value<bool>("loop_model:fiser_bb_loop_rmsd_corrected") == true){

		if (!ENV::is_set("pose:io:pdb:r")){
			cerr << "ERROR: reference pdb structure not set\n";
			return -1;
		}

		PRODART::POSE::pose_shared_ptr ref_pdb = PRODART::POSE::new_pose();
		ifstream ref_protein_file(PRODART::ENV::get_option_value<string>("pose:io:pdb:r").c_str(), ios::in);

		if (ref_protein_file.is_open()){
			ref_pdb->loadPdb(ref_protein_file);
			ref_pdb->set_label(PRODART::ENV::get_option_value<string>("pose:io:pdb:r"));
		}
		else {
			cerr << "ERROR: can't open file: " << PRODART::ENV::get_option_value<string>("pose:io:pdb:r")
				                                        		   << endl;
			return -1;
		}
		ref_protein_file.close();



		const int overwrite = PRODART::ENV::get_option_value<int>("loop_model:insert_overwrite");

		const_residue_shared_ptr res = in_pdb->get_residue(PRODART::ENV::get_option_value<string>("loop_model:insert_resnum"),
				PRODART::ENV::get_option_value<char>("loop_model:insert_chain"));


		if (res){

			const int resnum = res->get_internal_residue_index();

			const double rmsd = PROTOCOLS::LOOP::get_fiser_bb_loop_rmsd_corrected(ref_pdb, in_pdb, resnum, overwrite);

			cout << "LOOP_RMSD:\t" << rmsd << endl;

			return 0;
		}

		else {
			cerr << "ERROR: residue not found" << endl;
		}
	}
	else {
		cerr << "Nothing to do it's up to you" << endl;
	}



	return 0;

}



void print_phi_psi(pose_shared_ptr pose_, int start, int end, ostream& output  ){
	for (int i = start; i <= end; i++){

		const double phi = pose_->get_phi(i);
		const double psi = pose_->get_psi(i);
		const double omega = pose_->get_omega_to_prev(i);

		residue_shared_ptr res = pose_->get_residue(i);

		output << res->get_trimmed_pdb_residue_index() << "\t"
				//<< res->get_chain()->getChainID() << "\t"
				<< res->get_type().get_label3() << "\t"
				<< radians_to_degrees(phi) << "\t"
				<< radians_to_degrees(psi) << "\t"
				<< radians_to_degrees(omega) << "\t";
		if (i == start){
			output	<< "start";
		}
		else if (i == end){
			output	<< "end";
		}
		else {
			output	<< "middle";
		}
		output	<< endl;

	}
}



