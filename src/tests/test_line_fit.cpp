//
//  test_line_fit.cpp
//  pd2
//
//  Created by jmacdon on 02/01/2015.
//  Copyright (c) 2015 James. All rights reserved.
//

#include <stdio.h>
#include "backbonebuilder/alphabet_bb_builder.h"

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
#include "prodart_env/prodart_env.h"

#include "protocols/protocols.h"

#include "pose_utils/residue_reconstructor.h"
#include "rotamers/rotamer_db.h"
#include "rotamers/sidechain_builders/sidechain_builder.h"

#include "backbonebuilder/alphabet_bb_builder_manager.h"

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

using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace PRODART::ROTAMERS;
using namespace PRODART;
using namespace PRODART::POSE::BB_BUILDER;

int main( int argc, char *argv[] ) {
	PRODART::ENV::prodart_env::Instance()->init(argc, argv);

	if (argc != 2) {
		
		cerr << "Usage:\n"
		<< argv[0]<< " <PDB>\n"
		<< "Output:\n"
		<< "test"
		<< endl;
		
		return -1;
	}
	
	
	PRODART::POSE::pose_shared_ptr test_pose = PRODART::POSE::new_pose();
	ifstream protein_file(argv[1], ios::in);
	
	
	
	if (protein_file.is_open()){
		test_pose->loadPdb(protein_file);
		//validate_bb_pose(test_pose);
		validate_ca_pose(test_pose);
	}
	else {
		cerr << "ERROR: can't open file: " << argv[1]
		<< endl;
		return -1;
	}
	protein_file.close();
	
	PRODART::UTILS::vector3d_vector ca_coords;
	int size = test_pose->get_residue_count();
	
	for (int ii = 0; ii < size; ii++ ){
		ca_coords.push_back(test_pose->get_bb_atom(POSE::CA, ii)->get_coords());
	}
	PRODART::UTILS::vector3d first ,end;
	
	UTILS::line_fit3d(ca_coords, first, end);
	
	cout << "GSL_FIRST:\t" << first << endl;
	cout << "GSL_ENDL:\t" << end << endl << endl;
	
	PRODART::UTILS::vector3d eig3_first ,eig3_end;

	UTILS::eigen3_line_fit3d(ca_coords, eig3_first, eig3_end);

	cout << "EIGEN3_FIRST:\t" << eig3_first << endl;
	cout << "EIGEN3_ENDL:\t" << eig3_end << endl << endl;
	
}


