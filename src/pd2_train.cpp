/*
 * pd2_train.cpp
 *
 *  Created on: 21 Oct 2013
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
using boost::trim;

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

//! train potentials
int main( int argc, char *argv[] ) {

	PRODART::ENV::register_option("training_list",string(""),"input file list");
	PRODART::ENV::register_option("background_list",string(""),"input file list");
	PRODART::ENV::register_option("potential_to_train",string(""),"which potential to parameterise");
	PRODART::ENV::register_option("help",bool(false),"print options");

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


	if (ENV::is_set("training_list") && ENV::is_set("potential_to_train")){
		POTENTIALS::potential_shared_ptr pot = POTENTIALS::potentials_factory::Instance()->get_potential(POTENTIALS::potentials_name( PRODART::ENV::get_option_value<string>("potential_to_train") ));

		ifstream train_set(PRODART::ENV::get_option_value<string>("training_list").c_str(), ios::in);
	    string lineStr;
	    unsigned long length, lineNum = 0 ;
	    vector<string> SplitVec;
	    int pdbs_loaded = 0;

		cout << "reading file..." << endl;
		if (train_set.is_open()){
		    while ( !train_set.eof() ) {
		            getline(train_set, lineStr);

		            trim(lineStr);
		            lineNum++;

		            length = lineStr.length();

		            //cout << endl << lineNum << " " << length << " ";

		            if (length > 0) {

		                split( SplitVec, lineStr, is_any_of("\t ") );
		                string this_pdb = SplitVec[0];

		            	PRODART::POSE::pose_shared_ptr in_pose;

		            	bool loadOK = true;
		            	try {
		            		cout << "output_centred_backbone_segments: parsing: " << this_pdb << " ... ";
		            		in_pose = MISC::verified_load_bb_pose(this_pdb);//PRODART::POSE::new_pose();
		            	}
		            	catch (std::exception &e) {
		        			cout << "failed" << endl;
		            		std::cerr << "output_centred_backbone_segments: WARNING: skipping file: " << lineStr << endl;
		            		loadOK = false;
		            	}

		            	if (loadOK){
		            		ca_pose_meta_shared_ptr meta_ = new_ca_pose_meta(in_pose);
		            		if (pot->training_pose(meta_)){
		            			cout << "done" << endl;
		            			pdbs_loaded++;
		            		}
		            		else {
		            			cout << "failed" << endl;
		            		}
		            	}



		            }

		    }

		    cout << "pd2_train: INFO: " << pdbs_loaded << " PDB files have been successfully parsed" << endl;

		    pot->calc_potentials_using_training();
		    pot->output_training_results(cout);


		}
		else {
			cerr << "ERROR: can not open training set list" << ENV::get_option_value<string>("training_list") << endl;
		}
	}
	else {
		cerr << "ERROR: at least --training_list and --potential_to_train options need to be set" << endl;
	}




}


