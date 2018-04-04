/*
 * ca_bond.cpp
 *
 *  Created on: 26 Aug 2010
 *      Author: jmacdona
 */
#include "ca_bond.h"

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



potential_shared_ptr new_ca_bond(){
	potential_shared_ptr ptr(new ca_bond());
	return ptr;
}

ca_bond::ca_bond()  {

	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_bond"));

    trans_CA_CA_l0 = 3.8;
    cis_CA_CA_l0 = 2.95;

    trans_CA_CA_k = 890.0;
    cis_CA_CA_k = 85.0;


}

bool ca_bond::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:ca_bond");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();

	return true;
}

/*
bool ca_bond::provides(const potentials_name& query_name) const{
	if (query_name == name){
		return true;
	}

	return false;
}
*/


double ca_bond::bond_energy(PRODART::POSE::four_state_sec_struct secs, const double bond_length) const{

	double total_score = 0;

	if (secs == PRODART::POSE::ss4_CIS){
		total_score = (cis_CA_CA_k * pow(bond_length - cis_CA_CA_l0, 2));
	}
	else {
		total_score = (trans_CA_CA_k * pow(bond_length - trans_CA_CA_l0, 2));
	}

	return total_score;
}

void ca_bond::bond_grad_mag(PRODART::POSE::four_state_sec_struct secs, const double bond_length, double& energy, double& grad_mag) const{
	//energy = 0;
	//grad_mag = 0;

	if (secs == PRODART::POSE::ss4_CIS){
		const double diff = bond_length - cis_CA_CA_l0;
		grad_mag = (cis_CA_CA_k * 2.0 * (diff));
		energy = (cis_CA_CA_k * pow(diff, 2));
	}
	else {
		const double diff = bond_length - trans_CA_CA_l0;
		grad_mag = (trans_CA_CA_k * 2.0 * (diff));
		energy = (trans_CA_CA_k * pow(diff, 2));
	}

}

double ca_bond::bond_energy_grad(PRODART::POSE::four_state_sec_struct secs,
		PRODART::POSE::const_atom_shared_ptr atm1,
		PRODART::POSE::const_atom_shared_ptr atm2,
		UTILS::vector3d_vector& grad,
		const double weight) const{

	double energy = 0, grad_mag = 0;
	const double bond_len = (atm1->get_coords() - atm2->get_coords()).mod();
	bond_grad_mag(secs, bond_len, energy, grad_mag);
	grad_mag *= weight;
	const vector3d unit_grad_vec = (atm1->get_coords() - atm2->get_coords()) / bond_len;
	const vector3d atm1_grad = grad_mag * unit_grad_vec;
	const vector3d atm2_grad = -grad_mag * unit_grad_vec;

	const int index1 = atm1->get_seq_num();
	const int index2 = atm2->get_seq_num();
	grad[index1] += atm1_grad;
	grad[index2] += atm2_grad;


	return energy;
}

double ca_bond::frag_bond_energy(const PRODART::POSE::META::frag4_element& ele) const{

	double energy = 0;

	energy += bond_energy(ele.frag_ss_class, (ele.CA_atoms[1]->get_coords() - ele.CA_atoms[2]->get_coords()).mod());

	if (ele.frag_pos == PRODART::POSE::META::fragNTERM || ele.frag_pos == PRODART::POSE::META::fragNandCTERM){
		energy += bond_energy(PRODART::POSE::ss4_OTHER, (ele.CA_atoms[0]->get_coords() - ele.CA_atoms[1]->get_coords()).mod());
	}
	if (ele.frag_pos == PRODART::POSE::META::fragCTERM || ele.frag_pos == PRODART::POSE::META::fragNandCTERM){
		energy += bond_energy(PRODART::POSE::ss4_OTHER, (ele.CA_atoms[2]->get_coords() - ele.CA_atoms[3]->get_coords()).mod());
	}


	return energy;
}

double ca_bond::frag_bond_energy_grad(const PRODART::POSE::META::frag4_element& ele, vector3d_vector& grad, const double weight ) const{


	double energy = 0;

	energy += bond_energy_grad(ele.frag_ss_class, ele.CA_atoms[1], ele.CA_atoms[2], grad, weight);

	if (ele.frag_pos == PRODART::POSE::META::fragNTERM || ele.frag_pos == PRODART::POSE::META::fragNandCTERM){
		energy += bond_energy_grad(PRODART::POSE::ss4_OTHER, ele.CA_atoms[0], ele.CA_atoms[1], grad, weight);
	}
	if (ele.frag_pos == PRODART::POSE::META::fragCTERM || ele.frag_pos == PRODART::POSE::META::fragNandCTERM){
		energy += bond_energy_grad(PRODART::POSE::ss4_OTHER, ele.CA_atoms[2], ele.CA_atoms[3], grad, weight);
	}


	return energy;


}


double ca_bond::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);

	double total_score = 0;

	frag4_vector& fragments = ca_meta_dat->get_fragments();

	frag4_vector::const_iterator iter;
	for (iter = fragments.begin(); iter != fragments.end(); iter++){
		total_score += this->frag_bond_energy(*iter);
	}

	return energies_map.add_energy_component(this->name_vector[0], total_score);
}

// TODO
double ca_bond::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	const const_pose_shared_ptr _pose = _pose_meta->get_pose();
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	vector3d_vector& grad = _pose_meta->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);

	double total_score = 0;

	frag4_vector& fragments = ca_meta_dat->get_fragments();

	frag4_vector::const_iterator iter;
	for (iter = fragments.begin(); iter != fragments.end(); iter++){
		total_score += this->frag_bond_energy_grad(*iter, grad, weight);
	}

	return energies_map.add_energy_component(this->name_vector[0], total_score);
}


std::istream& ca_bond::load_data( std::istream& input ){
	string lineStr;


	long length, lineNum = 0 ;


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
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 2 ){

				string paraName = SplitVec[0];
				trim(paraName);

				if ( paraName.compare("trans_CA_CA_bond") == 0 ) {
					double l = lexical_cast<double>(SplitVec[1]);
					double k = lexical_cast<double>(SplitVec[2]);
					trans_CA_CA_l0 = l;
					trans_CA_CA_k = k;
				}
				else if ( paraName.compare("cis_CA_CA_bond") == 0 ){
					double l = lexical_cast<double>(SplitVec[1]);
					double k = lexical_cast<double>(SplitVec[2]);
					cis_CA_CA_l0 = l;
					cis_CA_CA_k = k;
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


}
}
}
}
