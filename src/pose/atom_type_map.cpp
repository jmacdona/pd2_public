/*
 * atom_type_map.cpp
 *
 *  Created on: 1 Aug 2011
 *      Author: jmacdona
 */
#include "atom_type_map.h"


#include "residue_type_map.h"

using namespace boost;
using namespace std;



namespace PRODART {
namespace POSE {

namespace {
const atom_type_map *instance = NULL;
boost::once_flag once_flag = BOOST_ONCE_INIT;
}

void atom_type_map::Init(){
	if (!instance){
		instance = new const atom_type_map;
	}
}


const atom_type_map* atom_type_map::Instance(){
	boost::call_once(&atom_type_map::Init, once_flag);
	return instance;
}

atom_type_map::atom_type_map(){

	// comment out later:
	//cerr << "atom_type_map: initialising..." << endl;

	assert(ENV::is_command_line_parsed());


	const string db_path = PRODART::ENV::get_option_value<string>("database:path:atom_type_map");
	std::ifstream input(db_path.c_str(), ios::in);

    string lineStr;


    long length, lineNum = 0 ;


    vector<string> SplitVec;
    while ( !input.eof() ) {
    	getline(input, lineStr);
    	string resStr;
    	lineNum++;

    	length = lineStr.length();

    	//cout << endl << lineNum << " " << length << " ";

    	if (length > 0) {
    		split( SplitVec, lineStr, is_any_of("\t") );
    		if ( SplitVec[0].substr(0,1).compare("#") != 0
    				&& SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 2 ){
    			trim(SplitVec[0]);
    			const string name = SplitVec[0];
    			const int ff_id = lexical_cast<int>(SplitVec[1]);

    			str_to_ff_id_map[name] = ff_id;
    			ff_id_to_str_map[ff_id] = name;



    		}

    	}

    }


	input.close();
}


atom_type_map::~atom_type_map(){
	delete instance;
	instance = 0;
}


}
}
