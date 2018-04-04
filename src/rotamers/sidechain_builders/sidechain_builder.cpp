/*
 * sidechain_builder.cpp
 *
 *  Created on: 12 Jan 2011
 *      Author: jmacdona
 */

#include "sidechain_builder.h"

using namespace boost;
using namespace std;
using namespace boost::filesystem;
using namespace PRODART::POSE_UTILS;
using namespace PRODART::POSE;

namespace PRODART {
namespace ROTAMERS {



namespace {
sidechain_builder *instance = NULL;
boost::once_flag sb_once_flag = BOOST_ONCE_INIT;
}


void sidechain_builder::Init(){
	if (!instance){
		instance = new sidechain_builder;
	}
	//return instance;
}


const sidechain_builder* sidechain_builder::Instance(){

	boost::call_once(&sidechain_builder::Init, sb_once_flag);

	return instance;
}

sidechain_builder::sidechain_builder(){
	this->init_resdefs();
}

sidechain_builder::~sidechain_builder(){
	delete instance;
	instance = 0;
}


void sidechain_builder::add_def(boost::filesystem::path def_path){

	residue_reconstructor_shared_ptr rc = new_residue_reconstructor();

	std::ifstream input(def_path.string().c_str(), ios::in);

	rc->load_residue_def(input);

	residue_type rt = rc->get_residue_type();

	//TODO add check for dupes

	this->res_rescontr_map[rt] = rc;



}

void sidechain_builder::init_resdefs(){
	cout << "sidechain_builder: loading residue definitions" << endl;

	const string defsDirPath = PRODART::ENV::get_option_value<string>("database:path:residue_defs");
	const path dir(defsDirPath);

	if (exists(dir)){
		directory_iterator end;
		for( directory_iterator iter(dir) ; iter != end ; ++iter )
			if ( is_directory( iter->status() ) )
			{

			}
			else {
				path resdef_path(iter->path());
				cout << "sidechain_builder: Opening definition:\t" << resdef_path << endl;
				add_def(resdef_path);
				//load_preset(wts_path);

			}
	}
}


POSE_UTILS::const_residue_reconstructor_shared_ptr sidechain_builder::get_reconstructor( POSE::residue_type type) const{

	return ((this->res_rescontr_map.find(type) != res_rescontr_map.end()) ? (res_rescontr_map.find(type)->second) : const_residue_reconstructor_shared_ptr());

}


}
}




