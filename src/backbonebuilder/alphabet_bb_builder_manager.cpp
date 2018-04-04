/*
 * alphabet_bb_builder_manager.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: jmacdona
 */
#include "alphabet_bb_builder_manager.h"


using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;

namespace PRODART {
namespace POSE {
namespace BB_BUILDER {



namespace {
alphabet_bb_builder_manager_shared_ptr instance;
boost::once_flag once_flag = BOOST_ONCE_INIT;
}

void alphabet_bb_builder_manager::Init(){
	if (!instance){
		instance = alphabet_bb_builder_manager_shared_ptr(new alphabet_bb_builder_manager());
	}
}

const_alphabet_bb_builder_manager_shared_ptr alphabet_bb_builder_manager::Instance(){
	/*
	if (!instance){
		instance = alphabet_bb_builder_manager_shared_ptr(new alphabet_bb_builder_manager());
	}
	*/
	boost::call_once(&alphabet_bb_builder_manager::Init, once_flag);
	return instance;
}



alphabet_bb_builder_manager::alphabet_bb_builder_manager() : default_even_len(6), default_odd_len(5), default_len(default_even_len){
	{ //L6
	const string l6_db_path = PRODART::ENV::get_option_value<string>("database:path:alphabet_bb_builder:L6");
	std::ifstream l6_inputdb(l6_db_path.c_str(), ios::in);
	cout << "alphabet_bb_builder_manager: loading db" << endl;
	alphabet_bb_builder_shared_ptr l6_ab = new_alphabet_bb_builder();
	l6_ab->load_alphabet_pose(l6_inputdb);
	alphabets[6] = l6_ab;
	}

	if (!this->get(default_len)){
		cerr << "alphabet_bb_builder_manager: WARNING: default_len alphabet not set in database path" << endl;
	}
	if (!this->get(default_even_len)){
		cerr << "alphabet_bb_builder_manager: WARNING: default_even_len alphabet not set in database path" << endl;
	}
	if (!this->get(default_odd_len)){
		cerr << "alphabet_bb_builder_manager: WARNING: default_odd_len alphabet not set in database path" << endl;
	}

}




alphabet_bb_builder_manager::~alphabet_bb_builder_manager(){

}



const_alphabet_bb_builder_shared_ptr alphabet_bb_builder_manager::get(const int alphabet_len) const{

	return alphabets.find(alphabet_len) != alphabets.end() ? (alphabets.find(alphabet_len))->second : const_alphabet_bb_builder_shared_ptr();

}

POSE::atom_shared_ptr_vector alphabet_bb_builder_manager::build_backbone( pose_shared_ptr pose_) const{
	return this->get(default_len)->build_backbone(pose_);
}

POSE::atom_shared_ptr_vector alphabet_bb_builder_manager::build_peptide_bonds(const int first_res, const int last_res, pose_shared_ptr pose_) const{
	return this->get(default_even_len)->build_peptide_bonds(first_res, last_res, pose_);
}












}
}
}




