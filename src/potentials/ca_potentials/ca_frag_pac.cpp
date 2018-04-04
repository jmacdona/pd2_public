/*
 * ca_frag_pac.cpp
 *
 *  Created on: Feb 27, 2011
 *      Author: jmacdon
 */

#include "ca_frag_pac.h"



using namespace boost;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;
using namespace std;

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

//const int NO_FRAG = 99999;

potential_shared_ptr new_ca_frag_pac(){
	potential_shared_ptr ptr(new ca_frag_pac());
	return ptr;
}



ca_frag_pac::ca_frag_pac()  {

	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_frag_pac"));

	int num = 2;
	restype_to_int_map[residue_type("ALA")] = num++;
	restype_to_int_map[residue_type("VAL")] = num++;
	restype_to_int_map[residue_type("PHE")] = num++;
	restype_to_int_map[residue_type("PRO")] = num++;
	restype_to_int_map[residue_type("MET")] = num++;
	restype_to_int_map[residue_type("ILE")] = num++;
	restype_to_int_map[residue_type("LEU")] = num++;
	restype_to_int_map[residue_type("ASP")] = num++;
	restype_to_int_map[residue_type("GLU")] = num++;
	restype_to_int_map[residue_type("LYS")] = num++;
	restype_to_int_map[residue_type("ARG")] = num++;
	restype_to_int_map[residue_type("SER")] = num++;
	restype_to_int_map[residue_type("THR")] = num++;
	restype_to_int_map[residue_type("TYR")] = num++;
	restype_to_int_map[residue_type("HIS")] = num++;
	restype_to_int_map[residue_type("CYS")] = num++;
	restype_to_int_map[residue_type("ASN")] = num++;
	restype_to_int_map[residue_type("GLN")] = num++;
	restype_to_int_map[residue_type("TRP")] = num++;
	restype_to_int_map[residue_type("GLY")] = num++;

}

bool ca_frag_pac::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:ca_frag_pac");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	return true;
}


std::istream& ca_frag_pac::load_data( std::istream& input ){

	string lineStr;

	long length, lineNum = 0 ;

	PAC_score_map.clear();

	string_vector SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 3 ){

				string paraName = SplitVec[0];
				trim(paraName);

				if ( paraName.compare("PAC_SS_SCORE") == 0 ){
					int PAC_bin = lexical_cast<int>(SplitVec[1]);
					double score = lexical_cast<double>(SplitVec[2]);
					PAC_score_map[PAC_bin] = score;
				}
				else {
                    cout << "ERROR - unknown parameter name: " << paraName << endl;
                    //fatalError = true;
                    //cout << paraValue << endl;
                }



			}
		}
	}


	return input;

}


int ca_frag_pac::rt_to_int(const POSE::residue_type &rt) const{
	return restype_to_int_map.find(rt) != restype_to_int_map.end() ? restype_to_int_map.find(rt)->second : 0;
}

int ca_frag_pac::overall_bin(const int frag_num,
			const POSE::residue_type &rt1,
			const POSE::residue_type &rt2) const{
	const int bin = frag_num
									+ 100 * rt_to_int(rt1)
									+ 10000 * rt_to_int(rt2);
	return bin;
}



double ca_frag_pac::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{

	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();

	double total_score = 0;

	//int frag_num;

	frag4_vector& fragments = ca_meta_dat->get_fragments();

	//frag4_vector::const_iterator iter;
	for (unsigned int i = 0; i < fragments.size(); i++){
		const int frag_num = fragments[i].frag_type_num;
		const int res_num1 = fragments[i].CA_1_residue_number;
		const residue_type rt1 = _pose->get_residue(res_num1)->get_type();
		const residue_type rt2 = _pose->get_residue(res_num1+1)->get_type();
		const int PAC_overall_profile = overall_bin(frag_num, rt1, rt2);

		const double this_energy = PAC_score_map.find(PAC_overall_profile) != PAC_score_map.end() ?  PAC_score_map.find(PAC_overall_profile)->second : 0.0;

		total_score += this_energy;

	}

    //cout << "BLAH!!!" << fragments.size() << "\t" << total_score << endl;

	return energies_map.add_energy_component(this->name_vector[0], total_score);
}

double ca_frag_pac::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	return get_energy(_pose_meta, energies_map);
}


}
}
}
}


