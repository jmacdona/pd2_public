/*
 * ca_frag_freq.cpp
 *
 *  Created on: 27 Aug 2010
 *      Author: jmacdona
 */

#include "ca_frag_freq.h"



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

const int NO_FRAG = 99999;

potential_shared_ptr new_ca_frag_freq(){
	potential_shared_ptr ptr(new ca_frag_freq());
	return ptr;
}



ca_frag_freq::ca_frag_freq()  {

	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_frag_freq"));

}

bool ca_frag_freq::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:ca_frag_freq");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	return true;
}

/*
bool ca_frag_freq::provides(const potentials_name& query_name) const{
	if (query_name == name){
		return true;
	}

	return false;
}
*/

std::istream& ca_frag_freq::load_data( std::istream& input ){

	string lineStr;

	long length, lineNum = 0 ;

	frag_energy.clear();
	frag12_energy.clear();
	frag13_energy.clear();
	frag123_energy.clear();

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

				if ( paraName.compare("frag_energy") == 0 ){

					const int frag_num = lexical_cast<int>(SplitVec[1]);

					const double this_frag_energy = lexical_cast<double>(SplitVec[2]);


					frag_energy[frag_num] = this_frag_energy;



				}
				else if ( paraName.compare("frag12_energy") == 0 ){
					const int frag_num = lexical_cast<int>(SplitVec[1]);

					const double this_frag_energy = lexical_cast<double>(SplitVec[2]);


					frag12_energy[frag_num] = this_frag_energy;
				}
				else if ( paraName.compare("frag13_energy") == 0 ){
					const int frag_num = lexical_cast<int>(SplitVec[1]);

					const double this_frag_energy = lexical_cast<double>(SplitVec[2]);


					frag13_energy[frag_num] = this_frag_energy;
				}
				else if (paraName.compare("frag123_energy") == 0 ){
					const int frag_num = lexical_cast<int>(SplitVec[1]);

					const double this_frag_energy = lexical_cast<double>(SplitVec[2]);


					frag123_energy[frag_num] = this_frag_energy;
				}
				else {
                    cout << "ERROR - unknown parameter name: " << paraName << endl;
                    //fatalError = true;
                    //cout << paraValue << endl;
                }


			}
		}
	}

	total_num_frags = frag_energy.size();
	return input;

}

int ca_frag_freq::get_combined_bin(const int n1, const int n2) const{
	return n1 + (this->total_num_frags * n2);
}

void ca_frag_freq::deconvolute_combined_bin(const int bin,int& n1, int& n2) const{
	n2 = bin / this->total_num_frags;

	n1 = bin % this->total_num_frags;
}

int ca_frag_freq::get_combined_bin(const int n1,
		const int n2,
		const int n3) const{

	return n1
		+ (this->total_num_frags * n2)
		+ (this->total_num_frags * this->total_num_frags * n3);
}

void ca_frag_freq::deconvolute_combined_bin(const int bin,
		int& n1,
		int& n2,
		int& n3) const{

	n3 = bin / (this->total_num_frags * this->total_num_frags);

	n2 = (bin - (n3 * this->total_num_frags * this->total_num_frags)) / this->total_num_frags;

	n1 = (bin - (n3 * this->total_num_frags * this->total_num_frags)) % this->total_num_frags;


}


double ca_frag_freq::get_frag_energy(const int prev_prev_frag_num, const int prev_frag_num, const int frag_num) const{

	double total_energy = 0;

    //const int frag_num = ele.frag_type_num;

	/*
    if (ele.frag_pos == fragNTERM || ele.frag_pos == fragNandCTERM){
    	prev_frag_num = NO_FRAG;
    	prev_prev_frag_num = NO_FRAG;
    }
    */
    if (prev_frag_num != NO_FRAG){
    	const int comb12_bin = this->get_combined_bin(prev_frag_num, frag_num);
    	total_energy += frag12_energy.find(comb12_bin) != frag12_energy.end() ? frag12_energy.find(comb12_bin)->second : 0.0;
    }
    if (prev_prev_frag_num != NO_FRAG){
    	const int comb13_bin = this->get_combined_bin(prev_prev_frag_num, frag_num);
    	total_energy += frag13_energy.find(comb13_bin) != frag13_energy.end() ? frag13_energy.find(comb13_bin)->second : 0.0; //this->frag13_energy[comb13_bin];
    	const int comb123_bin = this->get_combined_bin(prev_prev_frag_num, prev_frag_num, frag_num);
    	total_energy += frag123_energy.find(comb123_bin) != frag123_energy.end() ? frag123_energy.find(comb123_bin)->second : 0.0; //this->frag123_energy[comb123_bin];
    }

    total_energy += frag_energy.find(frag_num) != frag_energy.end() ? frag_energy.find(frag_num)->second : 0.0; //this->frag_energy[frag_num];

	return total_energy;
}


double ca_frag_freq::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{

	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();

	double total_score = 0;

	int prev_prev_frag_num = NO_FRAG;
	int prev_frag_num = NO_FRAG;
	int frag_num = NO_FRAG;

	frag4_vector& fragments = ca_meta_dat->get_fragments();

	//frag4_vector::const_iterator iter;
	for (unsigned int i = 0; i < fragments.size(); i++){
		frag_num = fragments[i].frag_type_num;
		if (fragments[i].frag_pos == fragNTERM || fragments[i].frag_pos == fragNandCTERM){
			prev_prev_frag_num = NO_FRAG;
			prev_frag_num = NO_FRAG;
		}

		total_score += this->get_frag_energy(prev_prev_frag_num, prev_frag_num, frag_num);

		prev_prev_frag_num = prev_frag_num;
		prev_frag_num = frag_num;

	}

    //cout << "BLAH!!!" << fragments.size() << "\t" << total_score << endl;

	return energies_map.add_energy_component(this->name_vector[0], total_score);
}

double ca_frag_freq::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	//const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	return get_energy(_pose_meta, energies_map);
}


}
}
}
}
