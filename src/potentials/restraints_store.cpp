/*
 * restraints_store.cpp
 *
 *  Created on: 25 Jan 2011
 *      Author: jmacdona
 */
#include "restraints_store.h"

#include "pose_meta/pose_meta_interface.h"

using namespace boost::filesystem;
using namespace std;

using std::string;
using std::ostream;
using std::istream;

using std::ofstream;

using std::ios;

using std::cout;
using std::endl;

using std::cerr;

using std::vector;
using std::cin;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::strncmp;
using std::exp;
using std::sqrt;
using boost::trim;
using boost::lexical_cast;
using boost::split;
using boost::is_any_of;

typedef vector<string> string_vector;


namespace PRODART {
namespace POSE {
namespace POTENTIALS {

namespace {
restraints_store *instance = NULL;
}


restraints_store::restraints_store(){

	if (ENV::is_set("restraints::rstfile")){
		cout << "restraints_store: loading restraints from command line..." << endl;
		const string db_path = ENV::get_option_value<string>("restraints::rstfile");
		std::ifstream inputdb(db_path.c_str(), ios::in);
		if (inputdb.is_open()){
			this->load_restraints(inputdb);
			cout << "restraints_store: ...restraints loaded\n" << endl;
		}
		else {
			cerr << "restraints_store: ERROR: can't open file: " << ENV::get_option_value<string>("restraints::rstfile") << endl;
			throw restraints_init_exception();
		}
		inputdb.close();

	}

}

restraints_store::~restraints_store(){
	delete instance;
	instance = NULL;
}

restraints_store* restraints_store::Instance(){
	if (!instance){
		instance = new restraints_store;
	}
	return instance;
}

void restraints_store::add_rst(const std::string& lineStr, const long lineNum){
	const long length = lineStr.length();
	string_vector SplitVec;
	string lineNumStr(" at line: ");
	if (lineNum >= 0){
		lineNumStr.append(lexical_cast<string>(lineNum));
	}
	else {
		lineNumStr.assign(": ");
	}
	if (length > 0) {
		split( SplitVec, lineStr, is_any_of("\t ") );
		if ( SplitVec[0].substr(0,1).compare("#") != 0
				&& SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 4 ){

			string paraName = SplitVec[0];

			if (paraName.compare("BOND_HARMONIC") == 0){
				if (SplitVec.size() < 8){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				bond_harmonic_rst_element ele;
				for (int i = 1; i <= 2; i++){
					const int start = 1 + ((i-1)*3);
					ele.resids[i-1] = SplitVec[start];
					ele.chains[i-1] = lexical_cast<char>(SplitVec[start+1]);
					ele.at_types[i-1] = POSE::atom_type(SplitVec[start+2]);
					cout << ele.resids[i-1] << "\t"
							<< ele.chains[i-1] << "\t"
							<< ele.at_types[i-1].get_label() << "\t";
				}
				const int eq_start = 1 + ((3-1)*3);
				ele.equil_val = lexical_cast<double>(SplitVec[eq_start]);
				ele.weight = lexical_cast<double>(SplitVec[eq_start+1]);
				cout << ele.equil_val << "\t"
						<< ele.weight << "\t";
				cout << endl;
				bond_rst_vec.push_back(ele);
			}
			else if (paraName.compare("ANGLE_HARMONIC") == 0){
				if (SplitVec.size() < 11){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				ang_harmonic_rst_element ele;
				for (int i = 1; i <= 3; i++){
					const int start = 1 + ((i-1)*3);
					ele.resids[i-1] = SplitVec[start];
					ele.chains[i-1] = lexical_cast<char>(SplitVec[start+1]);
					ele.at_types[i-1] = POSE::atom_type(SplitVec[start+2]);
					cout << ele.resids[i-1] << "\t"
							<< ele.chains[i-1] << "\t"
							<< ele.at_types[i-1].get_label() << "\t";
				}
				const int eq_start = 1 + ((4-1)*3);
				ele.equil_val = UTILS::degrees_to_radians(lexical_cast<double>(SplitVec[eq_start]));
				ele.weight = lexical_cast<double>(SplitVec[eq_start+1]);
				cout << ele.equil_val << "\t"
						<< ele.weight << "\t";
				cout << endl;
				ang_rst_vec.push_back(ele);
			}
			else if (paraName.compare("DIHEDRAL_HARMONIC") == 0){
				if (SplitVec.size() < 14){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				dih_harmonic_rst_element ele;
				for (int i = 1; i <= 4; i++){
					const int start = 1 + ((i-1)*3);
					ele.resids[i-1] = SplitVec[start];
					ele.chains[i-1] = lexical_cast<char>(SplitVec[start+1]);
					ele.at_types[i-1] = POSE::atom_type(SplitVec[start+2]);
					cout << ele.resids[i-1] << "\t"
							<< ele.chains[i-1] << "\t"
							<< ele.at_types[i-1].get_label() << "\t";
				}
				const int eq_start = 1 + ((5-1)*3);
				ele.equil_val = UTILS::degrees_to_radians(lexical_cast<double>(SplitVec[eq_start]));
				ele.weight = lexical_cast<double>(SplitVec[eq_start+1]);
				cout << ele.equil_val << "\t"
						<< ele.weight << "\t";
				cout << endl;
				dih_rst_vec.push_back(ele);
			}
			else if (paraName.compare("SEC_STRUCT") == 0){
				if (SplitVec.size() < 5){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				sec_struct_rst_element ele;


				ele.resids[0] = SplitVec[1];
				ele.chains[0] = lexical_cast<char>(SplitVec[2]);

				cout << ele.resids[0] << "\t"
						<< ele.chains[0] << "\t";

				std::string secstring = SplitVec[3];

				if (secstring.compare("H") == 0){
					ele.secs = POSE::ss4_HELIX;
				}
				else if (secstring.compare("E") == 0){
					ele.secs = POSE::ss4_STRAND;
				}
				else if (secstring.compare("-") == 0){
					ele.secs = POSE::ss4_OTHER;
				}
				else if (secstring.compare("C") == 0){
					ele.secs = POSE::ss4_CIS;
				}
				else {
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}

				ele.weight = lexical_cast<double>(SplitVec[4]);
				cout << POSE::to_char(ele.secs) << "\t"
						<< ele.weight;
				cout << endl;
				secs_rst_vec.push_back(ele);

			}
			else if (paraName.compare("STRAND_PAIR") == 0){
				if (SplitVec.size() < 7){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				bond_harmonic_rst_element ele;


				ele.resids[0] = SplitVec[1];
				ele.chains[0] = lexical_cast<char>(SplitVec[2]);
				ele.resids[1] = SplitVec[3];
				ele.chains[1] = lexical_cast<char>(SplitVec[4]);

				cout << ele.resids[0] << "\t"
						<< ele.chains[0] << "\t";
				cout << ele.resids[1] << "\t"
						<< ele.chains[1] << "\t";

				const int eq_start = 5;
				if (SplitVec[eq_start].compare("PARALLEL") == 0){
					ele.is_parallel = true;
				}
				else if (SplitVec[eq_start].compare("ANTIPARALLEL") == 0){
					ele.is_parallel = false;
				}
				else {
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				ele.weight = lexical_cast<double>(SplitVec[eq_start+1]);
				cout << SplitVec[eq_start] << "\t"
						<< ele.weight << "\t";
				cout << endl;
				strand_pair_vec.push_back(ele);
			}
			else if (paraName.compare("COORD_HARMONIC") == 0){
				if (SplitVec.size() < 9){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				coord_rst_element ele;
				for (int i = 1; i <= 1; i++){
					const int start = 1 + ((i-1)*3);
					ele.resids[i-1] = SplitVec[start];
					ele.chains[i-1] = lexical_cast<char>(SplitVec[start+1]);
					ele.at_types[i-1] = POSE::atom_type(SplitVec[start+2]);
					cout << ele.resids[i-1] << "\t"
							<< ele.chains[i-1] << "\t"
							<< ele.at_types[i-1].get_label() << "\t";
				}
				const int eq_start = 1 + ((2-1)*3);
				const double x = lexical_cast<double>(SplitVec[eq_start]);
				const double y = lexical_cast<double>(SplitVec[eq_start+1]);
				const double z = lexical_cast<double>(SplitVec[eq_start+2]);
				ele.coord = UTILS::vector3d(x,y,z);
				ele.weight = lexical_cast<double>(SplitVec[eq_start+3]);
				ele.equil_val = lexical_cast<double>(SplitVec[eq_start+4]);

				cout << ele.coord.x << "\t"
						 << ele.coord.y << "\t"
						 << ele.coord.z << "\t"
						<< ele.weight << "\t"
						<< ele.equil_val << "\t";
				cout << endl;
				coord_rst_vec.push_back(ele);
			}
			else if (paraName.compare("RESIDUE_CONTACT") == 0){
				if (SplitVec.size() < 8){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				contact_rst_element ele;


				ele.resids[0] = SplitVec[1];
				ele.chains[0] = lexical_cast<char>(SplitVec[2]);
				ele.resids[1] = SplitVec[3];
				ele.chains[1] = lexical_cast<char>(SplitVec[4]);

				cout << ele.resids[0] << "\t"
						<< ele.chains[0] << "\t";
				cout << ele.resids[1] << "\t"
						<< ele.chains[1] << "\t";

				const int eq_start = 5;

				ele.equil_val_lower = lexical_cast<double>(SplitVec[eq_start]);
				ele.equil_val_upper = lexical_cast<double>(SplitVec[eq_start+1]);
				ele.weight = lexical_cast<double>(SplitVec[eq_start+2]);


				cout << ele.equil_val_lower << "\t"
						<< ele.equil_val_upper << "\t"
						<< ele.weight << "\t";
				cout << endl;

				if (ele.equil_val_upper < ele.equil_val_lower){
					cerr << "restraints_store: load_restraints: ERROR: RESIDUE_CONTACT restraint lower bound: " <<  ele.equil_val_lower
							<< " must be lower than upper bound: " << ele.equil_val_upper
							<< endl;
					throw restraints_init_exception();
				}

				contact_rst_vec.push_back(ele);
			}
			else if (paraName.compare("CA_GO_CONTACT") == 0){
				if (SplitVec.size() < 8){
					cerr << "restraints_store: load_restraints: rstfile format error" << lineNumStr << " parameter: " << paraName << endl;
					throw restraints_init_exception();
				}
				cout << paraName << "\t";
				contact_rst_element ele;


				ele.resids[0] = SplitVec[1];
				ele.chains[0] = lexical_cast<char>(SplitVec[2]);
				ele.resids[1] = SplitVec[3];
				ele.chains[1] = lexical_cast<char>(SplitVec[4]);

				cout << ele.resids[0] << "\t"
						<< ele.chains[0] << "\t";
				cout << ele.resids[1] << "\t"
						<< ele.chains[1] << "\t";

				const int eq_start = 5;

				ele.equil_val_lower = lexical_cast<double>(SplitVec[eq_start]);
				ele.equil_val_upper = lexical_cast<double>(SplitVec[eq_start+1]);
				ele.weight = lexical_cast<double>(SplitVec[eq_start+2]);


				cout << ele.equil_val_lower << "\t"
						<< ele.equil_val_upper << "\t"
						<< ele.weight << "\t";
				cout << endl;

				if (ele.equil_val_upper < ele.equil_val_lower){
					cerr << "restraints_store: load_restraints: ERROR: RESIDUE_CONTACT restraint lower bound: " <<  ele.equil_val_lower
							<< " must be lower than upper bound: " << ele.equil_val_upper
							<< endl;
					throw restraints_init_exception();
				}


				ca_GO_rst_vec.push_back(ele);
			}
			else {
				cerr << "restraints_store: load_restraints: ERROR - unknown parameter name" << lineNumStr << " parameter: " << paraName << endl;
				throw restraints_init_exception();
			}

		}
	}
}

std::istream& restraints_store::load_restraints(std::istream& input){

	string lineStr;


	//long length,
	long lineNum = 0 ;


	string_vector SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		//length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";
		add_rst(lineStr, lineNum);

	}


	return input;
}

void restraints_store::apply_restraints(PRODART::POSE::META::pose_meta_shared_ptr pose_meta_){
	if (bond_rst_vec.size() != 0
			|| ang_rst_vec.size() != 0
			|| dih_rst_vec.size() != 0
			|| strand_pair_vec.size() != 0
			|| sse_axis_rst_vec.size() != 0
			|| ca_GO_rst_vec.size() != 0
			|| coord_rst_vec.size() != 0 ){
		cout << "restraints_store::apply_restraints: applying restraints..." << endl;

		POSE::pose_shared_ptr pose_ = pose_meta_->get_pose();

		unsigned long rsts_added = 0;

		for (bond_harm_rst_vector::const_iterator it = bond_rst_vec.begin(); it != bond_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			vector<atom_shared_ptr> atoms(it->get_n());
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}
				atoms[i] = res[i]->get_atom(it->at_types[i]);
				if (!atoms[i]){
					cerr << "restraints_store::apply_restraints: Warning: atom not found in residue: " << it->at_types[i].get_label() << "\t"
							<< it->chains[i] << endl;
				}
			}
			rsts_added++;
			pose_meta_->add_dist_harmonic_restraint(atoms[0]->get_seq_num(), atoms[1]->get_seq_num(), it->equil_val, it->weight);

		}

		for (sse_axis_rst_vector::const_iterator it = sse_axis_rst_vec.begin(); it != sse_axis_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			vector<atom_shared_ptr> atoms(it->get_n());
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}
				atoms[i] = res[i]->get_atom(it->at_types[i]);
				if (!atoms[i]){
					cerr << "restraints_store::apply_restraints: Warning: atom not found in residue: " << it->at_types[i].get_label() << "\t"
							<< it->chains[i] << endl;
				}
			}
			rsts_added++;
			pose_meta_->add_sse_restraint(res[0]->get_internal_residue_index(), res[1]->get_internal_residue_index(), it->secs, it->weight, it->coord, it->coord2);
			//pose_meta_->add_harmonic_restraint(atoms[0]->get_seq_num(), atoms[1]->get_seq_num(), it->equil_val, it->weight);

		}

		for (ang_harm_rst_vector::const_iterator it = ang_rst_vec.begin(); it != ang_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			vector<atom_shared_ptr> atoms(it->get_n());
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}
				atoms[i] = res[i]->get_atom(it->at_types[i]);
				if (!atoms[i]){
					cerr << "restraints_store::apply_restraints: Warning: atom not found in residue: " << it->at_types[i].get_label() << "\t"
							<< it->chains[i] << endl;
				}
			}
			rsts_added++;
			pose_meta_->add_angle_harmonic_restraint(atoms[0]->get_seq_num(),
					atoms[1]->get_seq_num(),
					atoms[2]->get_seq_num(),
					it->equil_val, it->weight);

		}

		for (dih_harm_rst_vector::const_iterator it = dih_rst_vec.begin(); it != dih_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			vector<atom_shared_ptr> atoms(it->get_n());
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}
				atoms[i] = res[i]->get_atom(it->at_types[i]);
				if (!atoms[i]){
					cerr << "restraints_store::apply_restraints: Warning: atom not found in residue: " << it->at_types[i].get_label() << "\t"
							<< it->chains[i] << endl;
				}
			}
			rsts_added++;
			pose_meta_->add_dihedral_harmonic_restraint(atoms[0]->get_seq_num(),
					atoms[1]->get_seq_num(),
					atoms[2]->get_seq_num(),
					atoms[3]->get_seq_num(),
					it->equil_val, it->weight);

		}

		for (sec_struct_rst_vector::const_iterator it = secs_rst_vec.begin(); it != secs_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			//atom_shared_ptr atoms[it->get_n()];
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}

			}
			rsts_added++;
			pose_meta_->add_sec_struct_restraint(res[0]->get_internal_residue_index(),
					it->secs, it->weight);

		}
		for (bond_harm_rst_vector::const_iterator it = strand_pair_vec.begin(); it != strand_pair_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			//atom_shared_ptr atoms[it->get_n()];
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}

			}
			rsts_added++;
			pose_meta_->add_strand_residue_pair_restraint(res[0]->get_internal_residue_index(), res[1]->get_internal_residue_index(), it->is_parallel, it->weight);

		}
		for (coord_rst_vector::const_iterator it = coord_rst_vec.begin(); it != coord_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			vector<atom_shared_ptr> atoms(it->get_n());
			//atom_shared_ptr atoms[it->get_n()];
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}
				atoms[i] = res[i]->get_atom(it->at_types[i]);
				if (!atoms[i]){
					cerr << "restraints_store::apply_restraints: Warning: atom not found in residue: " << it->at_types[i].get_label() << "\t"
							<< it->chains[i] << endl;
				}

			}
			rsts_added++;
			pose_meta_->add_harmonic_coord_restraint(atoms[0]->get_seq_num(),
					it->coord,
					it->weight,
					it->equil_val);

		}
		for (contact_rst_vector::const_iterator it = contact_rst_vec.begin(); it != contact_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			//atom_shared_ptr atoms[it->get_n()];
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}

			}


			rsts_added++;
			pose_meta_->add_residue_contact_restraint(res[0]->get_internal_residue_index(), res[1]->get_internal_residue_index(),
					it->equil_val_lower, it->equil_val_upper, it->weight);


		}

		for (contact_rst_vector::const_iterator it = ca_GO_rst_vec.begin(); it != ca_GO_rst_vec.end(); it++){
			vector<residue_shared_ptr> res(it->get_n());
			//atom_shared_ptr atoms[it->get_n()];
			for (int i = 0; i < it->get_n(); i++){
				res[i] = pose_->get_residue(it->resids[i], it->chains[i]);
				if (!res[i]){
					cerr << "restraints_store::apply_restraints: ERROR: residue does not exist: " << it->resids[i] << "\t"
							<< it->chains[i] << endl;
					throw restraints_init_exception();

				}

			}


			// TODO ADD TO META DAT

			rsts_added++;
			pose_meta_->add_ca_GO_restraint(res[0]->get_internal_residue_index(), res[1]->get_internal_residue_index(),
					it->equil_val_lower, it->equil_val_upper, it->weight);
			//cout << "restraints_store::apply_restraints: Added Go contact restraint." << endl;


		}


		cout << "restraints_store::apply_restraints: added " << rsts_added << " restraints" << endl;

	}
}


void restraints_store::add_sec_struct_rst(const std::string res1, const char ch1,
		const POSE::four_state_sec_struct secs, const double weight){
	sec_struct_rst_element ele;


	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.secs = secs;
	ele.weight = weight;

	secs_rst_vec.push_back(ele);
}
void restraints_store::add_coord_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
		const UTILS::vector3d pos, const double weight, const double equil_dist){
	coord_rst_element ele;
	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.at_types[0] = at1;
	ele.weight = weight;
	ele.equil_val = equil_dist;
	ele.coord = pos;

	coord_rst_vec.push_back(ele);
}
void restraints_store::add_bond_harm_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
		const std::string res2, const char ch2, const POSE::atom_type at2,
		const double equil_val, const double weight){
	bond_harmonic_rst_element ele;

	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.at_types[0] = at1;
	ele.resids[1] = res2;
	ele.chains[1] = ch2;
	ele.at_types[1] = at2;

	ele.equil_val = equil_val;
	ele.weight = weight;

	bond_rst_vec.push_back(ele);
}
void restraints_store::add_ang_harm_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
		const std::string res2, const char ch2, const POSE::atom_type at2,
		const std::string res3, const char ch3, const POSE::atom_type at3,
		const double equil_val, const double weight){
	ang_harmonic_rst_element ele;

	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.at_types[0] = at1;

	ele.resids[1] = res2;
	ele.chains[1] = ch2;
	ele.at_types[1] = at2;

	ele.resids[2] = res3;
	ele.chains[2] = ch3;
	ele.at_types[2] = at3;

	ele.equil_val = equil_val;
	ele.weight = weight;

	ang_rst_vec.push_back(ele);
}
void restraints_store::add_dih_harm_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
		const std::string res2, const char ch2, const POSE::atom_type at2,
		const std::string res3, const char ch3, const POSE::atom_type at3,
		const std::string res4, const char ch4, const POSE::atom_type at4,
		const double equil_val, const double weight){

	dih_harmonic_rst_element ele;

	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.at_types[0] = at1;

	ele.resids[1] = res2;
	ele.chains[1] = ch2;
	ele.at_types[1] = at2;

	ele.resids[2] = res3;
	ele.chains[2] = ch3;
	ele.at_types[2] = at3;

	ele.resids[3] = res4;
	ele.chains[3] = ch4;
	ele.at_types[3] = at4;

	ele.equil_val = equil_val;
	ele.weight = weight;

	dih_rst_vec.push_back(ele);
}
void restraints_store::add_strand_pair_rst(const std::string res1, const char ch1,
		const std::string res2, const char ch2,
		const bool is_parallel, const double weight){
	bond_harmonic_rst_element ele;


	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.resids[1] = res2;
	ele.chains[1] = ch2;



	ele.is_parallel = is_parallel;
	ele.weight = weight;

	strand_pair_vec.push_back(ele);
}

void restraints_store::add_sse_axis_rst(const std::string res1, const char ch1,
		const std::string res2, const char ch2,
		POSE::four_state_sec_struct sse_type, const double weight,
		const UTILS::vector3d start, const UTILS::vector3d end){
	sse_axis_rst_element ele;

	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.at_types[0] = POSE::atom_type("CA");
	ele.resids[1] = res2;
	ele.chains[1] = ch2;
	ele.at_types[1] = POSE::atom_type("CA");

	ele.secs = sse_type;
	ele.weight = weight;

	ele.coord = start;
	ele.coord2 = end;

	sse_axis_rst_vec.push_back(ele);
}

void restraints_store::add_contact_rst(const std::string res1, const char ch1,
		const std::string res2, const char ch2,
		const double equil_val_lower, const double equil_val_higher,
		const double weight){
	contact_rst_element ele;

	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.resids[1] = res2;
	ele.chains[1] = ch2;

	ele.equil_val_lower = equil_val_lower;
	ele.equil_val_upper = equil_val_higher;
	ele.weight = weight;

	if (ele.equil_val_upper < ele.equil_val_lower){
		cerr << "restraints_store: load_restraints: ERROR: RESIDUE_CONTACT restraint lower bound: " <<  ele.equil_val_lower
				<< " must be lower than upper bound: " << ele.equil_val_upper
				<< endl;
	}
	else {
		contact_rst_vec.push_back(ele);
	}
}

void restraints_store::add_ca_GO_rst(const std::string res1, const char ch1,
		const std::string res2, const char ch2,
		const double equil_val_lower, const double equil_val_higher,
		const double weight){
	contact_rst_element ele;

	ele.resids[0] = res1;
	ele.chains[0] = ch1;
	ele.resids[1] = res2;
	ele.chains[1] = ch2;

	ele.equil_val_lower = equil_val_lower;
	ele.equil_val_upper = equil_val_higher;
	ele.weight = weight;

	if (ele.equil_val_upper < ele.equil_val_lower){
		cerr << "restraints_store: load_restraints: ERROR: RESIDUE_CONTACT restraint lower bound: " <<  ele.equil_val_lower
				<< " must be lower than upper bound: " << ele.equil_val_upper
				<< endl;
	}
	else {
		ca_GO_rst_vec.push_back(ele);
	}
}

}
}
}
