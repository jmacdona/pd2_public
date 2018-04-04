/*
 * ca_tau_omega.cpp
 *
 *  Created on: 26 Aug 2010
 *      Author: jmacdona
 */
#include "ca_tau_omega.h"



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

potential_shared_ptr new_ca_tau_omega(){
	potential_shared_ptr ptr(new ca_tau_omega());
	return ptr;
}


ca_tau_omega::ca_tau_omega()  {

	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_tau_omega"));

}

bool ca_tau_omega::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:ca_tau_omega");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	return true;
}

/*
bool ca_tau_omega::provides(const potentials_name& query_name) const{
	if (query_name == name){
		return true;
	}

	return false;
}
*/

std::istream& ca_tau_omega::load_data( std::istream& input ){
	string lineStr;

	long length, lineNum = 0 ;

	tau_0_params.clear();
	tau_k_params.clear();
	tau_const_params.clear();
	omega1_0_params.clear();
	omega1_k_params.clear();
	omega1_const_params.clear();
	omega2_0_params.clear();
	omega2_k_params.clear();
	omega2_const_params.clear();

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
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 9 ){

				string paraName = SplitVec[0];
				trim(paraName);

				if ( paraName.compare("frag_params") == 0 ){

					const int frag_num = lexical_cast<int>(SplitVec[1]);

					const double tau_0 = lexical_cast<double>(SplitVec[3]);
					const double tau_k = lexical_cast<double>(SplitVec[4]);

					const double omega1_0 = lexical_cast<double>(SplitVec[5]);
					const double omega1_k = lexical_cast<double>(SplitVec[6]);

					const double omega2_0 = lexical_cast<double>(SplitVec[7]);
					const double omega2_k = lexical_cast<double>(SplitVec[8]);

					tau_0_params[frag_num] = tau_0;
					tau_k_params[frag_num] = tau_k;

					omega1_0_params[frag_num] = omega1_0;
					omega1_k_params[frag_num] = omega1_k;

					omega2_0_params[frag_num] = omega2_0;
					omega2_k_params[frag_num] = omega2_k;



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

double ca_tau_omega::get_dihedral_distance(const double dih1, const double dih2)const{

	const double diff = dih2 - dih1;

	if (fabs(diff) < PI){
		return diff;
	}
	else if (diff > 0){
		return diff - (2.0 * PI);
	}
	else {
		return diff + (2.0 * PI);
	}

}


double ca_tau_omega::frag_energy(const PRODART::POSE::META::frag4_element& ele) const{

	double total_energy = 0;

    const int frag_num = ele.frag_type_num;

    const double tau = ele.tau;
    const double omega1 = ele.omega1;
    const double omega2 = ele.omega2;

    const double tau_0 = tau_0_params.find(frag_num) != tau_0_params.end() ? tau_0_params.find(frag_num)->second : 0.0  ;//this->tau_0_params[frag_num];
    const double omega1_0 = omega1_0_params.find(frag_num) != omega1_0_params.end() ? omega1_0_params.find(frag_num)->second : 0.0  ;//this->omega1_0_params[frag_num];
    const double omega2_0 = omega2_0_params.find(frag_num) != omega2_0_params.end() ? omega2_0_params.find(frag_num)->second : 0.0  ;//this->omega2_0_params[frag_num];

    const double tau_k = tau_k_params.find(frag_num) != tau_k_params.end() ? tau_k_params.find(frag_num)->second : 0.0  ;//this->tau_k_params[frag_num];
    const double omega1_k = omega1_k_params.find(frag_num) != omega1_k_params.end() ? omega1_k_params.find(frag_num)->second : 0.0  ;//this->omega1_k_params[frag_num];
    const double omega2_k = omega2_k_params.find(frag_num) != omega2_k_params.end() ? omega2_k_params.find(frag_num)->second : 0.0  ;//this->omega2_k_params[frag_num];

    const double tau_score = tau_k * pow(this->get_dihedral_distance(tau, tau_0),2);

    const double omega1_score = omega1_k* pow(omega1 - omega1_0,2);

    const double omega2_score = omega2_k* pow(omega2 - omega2_0,2);


    if (ele.frag_pos == fragMIDDLE){
    	total_energy += 0.5 * (omega1_score + omega2_score);
    }
    else if (ele.frag_pos == fragNTERM){
    	total_energy += omega1_score + (0.5 * omega2_score);
    }
    else if (ele.frag_pos == fragCTERM){
    	total_energy += omega2_score + (0.5 * omega1_score);
    }
    else if (ele.frag_pos == fragNandCTERM){
    	total_energy += (0.5 * omega1_score) + (0.5 * omega2_score);
    }
    total_energy += tau_score;

    /*
    cout << frag_num << "\t"
    		<< ele.tau << "\t"
    		<<  ele.omega1 << "\t"
    		<< ele.omega2 << "\t"
    		<< tau_0 << "\t"
    		;
    cout << total_energy << endl;
    */

	return total_energy;
}


double ca_tau_omega::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{

	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();

	double total_score = 0;

	frag4_vector& fragments = ca_meta_dat->get_fragments();

	frag4_vector::const_iterator iter;
	for (iter = fragments.begin(); iter != fragments.end(); iter++){
		total_score += this->frag_energy(*iter);
	}

    //cout << "BLAH!!!" << fragments.size() << "\t" << total_score << endl;

	return energies_map.add_energy_component(this->name_vector[0], total_score);
}

double ca_tau_omega::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	//const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	return get_energy(_pose_meta, energies_map);
}



}
}
}
}

