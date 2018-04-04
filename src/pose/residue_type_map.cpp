//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * residue_type_map.cpp
 *
 *  Created on: Feb 20, 2010
 *      Author: jmacdon
 */

#include "residue_type_map.h"

using namespace boost;
using namespace std;



namespace PRODART {
namespace POSE {

namespace {
const residue_type_map *instance = NULL;
boost::once_flag once_flag = BOOST_ONCE_INIT;
}


void residue_type_map::Init(){
	if (!instance){
		instance = new const residue_type_map;
	}
}

const residue_type_map* residue_type_map::Instance(){
	/*
	if (!instance){
		instance = new const residue_type_map;
	}
	*/
	boost::call_once(&residue_type_map::Init, once_flag);
	return instance;
}

//const residue_type_map *residue_type_map::instance = NULL;

/*
enum StandardResidueType {UNDEF_RES, UNK_RES, ALA, VAL, PHE, PRO, MET, ILE, LEU, ASP, GLU, LYS, ARG, SER,
	THR, TYR, HIS, CYS, ASN, GLN, TRP, GLY};
*/

residue_type_map::residue_type_map(){

	assert(ENV::is_command_line_parsed());


	one_letter_name_map.clear();
	char_to_res_type_map.clear();

	one_letter_name_map[residue_type("ALA")] = 'A';
	one_letter_name_map[residue_type("VAL")] = 'V';
	one_letter_name_map[residue_type("PHE")] = 'F';
	one_letter_name_map[residue_type("PRO")] = 'P';
	one_letter_name_map[residue_type("MET")] = 'M';
	one_letter_name_map[residue_type("ILE")] = 'I';
	one_letter_name_map[residue_type("LEU")] = 'L';
	one_letter_name_map[residue_type("ASP")] = 'D';
	one_letter_name_map[residue_type("GLU")] = 'E';
	one_letter_name_map[residue_type("LYS")] = 'K';
	one_letter_name_map[residue_type("ARG")] = 'R';
	one_letter_name_map[residue_type("SER")] = 'S';
	one_letter_name_map[residue_type("THR")] = 'T';
	one_letter_name_map[residue_type("TYR")] = 'Y';
	one_letter_name_map[residue_type("HIS")] = 'H';
	one_letter_name_map[residue_type("CYS")] = 'C';
	one_letter_name_map[residue_type("ASN")] = 'N';
	one_letter_name_map[residue_type("GLN")] = 'Q';
	one_letter_name_map[residue_type("TRP")] = 'W';
	one_letter_name_map[residue_type("GLY")] = 'G';

	one_letter_name_map[residue_type("ASX")] = 'B';

	one_letter_name_map[residue_type("UNK")] = 'U';
	one_letter_name_map[residue_type()] = '-';

	residue_type_char_map::const_iterator iter;

	for (iter = one_letter_name_map.begin() ; iter != one_letter_name_map.end(); iter++){
		char_to_res_type_map[iter->second] = iter->first;
	}



}


residue_type_map::~residue_type_map(){
	delete instance;
	instance = 0;
}






residue_type_vector one_letter_seq_to_residue_type_vector(std::string seq){
	const residue_type_map  * rt_map = residue_type_map::Instance();
	trim(seq);
	residue_type_vector rtn_vec;
	rtn_vec.reserve(seq.size());
	for (unsigned int i = 0; i < seq.size(); i++){
		rtn_vec.push_back(rt_map->to_residue_type(seq.at(i)));
	}

	return rtn_vec;

}



}
}
