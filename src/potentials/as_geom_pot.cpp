/*
 * as_geom_pot.cpp
 *
 *  Created on: Jan 2, 2011
 *      Author: jmacdon
 */

#include "as_geom_pot.h"

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

potential_shared_ptr new_as_geom_pot(){
	potential_shared_ptr ptr(new as_geom_pot());
	return ptr;
}

//const prot_backbone_map* const as_geom_pot::bb_map(prot_backbone_map::Instance());
const double as_geom_pot::grad_h = 1e-8;


as_geom_pot::as_geom_pot() : bb_map(prot_backbone_map::Instance()){
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("as_geom_pot"));
	motif_loaded = false;
	motif_pdb = new_pose();
	using_fixed_mapping = false;

	// TODO: not ideal needs fixing
	rand_num = PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id());
}


void as_geom_pot::get_ca_atom_mapping(POSE::pose_shared_ptr pose_,
		std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> &atom_mapping,
		pose_shared_ptr this_motif) const{

	if (!this_motif){
		this_motif = motif_pdb;
	}

	atom_mapping.clear();

	for (unsigned int ele = 0; ele < motif_target_res_map_vector.size(); ele++){
		int_int_map::const_iterator iter;
		for (iter = motif_target_res_map_vector[ele].begin(); iter != motif_target_res_map_vector[ele].end(); iter++ ){

			atom_shared_ptr mot_at = this_motif->get_bb_atom(POSE::CA, iter->first);
			if (mot_at->isSet() && mot_at->isActive()){
				atom_shared_ptr pose_at = pose_->get_bb_atom(POSE::CA, iter->second);
				if (pose_at->isActive() && pose_at->isSet()){
					atom_mapping[mot_at] = pose_at;
				}
			}

		}
	}

}

double as_geom_pot::get_ca_rmsd(POSE::pose_shared_ptr pose_) const{
	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;
	this->get_ca_atom_mapping(pose_, atom_mapping);
	return POSE_UTILS::get_rmsd(atom_mapping);
}

void as_geom_pot::get_bb_atom_mapping(POSE::pose_shared_ptr pose_,
		std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> &atom_mapping,
		pose_shared_ptr this_motif) const{

	if (!this_motif){
		this_motif = motif_pdb;
	}

	const std::vector<BBAtomType> bb_atom_list = bb_map->get_small_bb_atom_list();

	atom_mapping.clear();


	for (unsigned int ele = 0; ele < motif_target_res_map_vector.size(); ele++){
		int_int_map::const_iterator iter;
		for (iter = motif_target_res_map_vector[ele].begin(); iter != motif_target_res_map_vector[ele].end(); iter++ ){
			for (std::vector<BBAtomType>::const_iterator at_it = bb_atom_list.begin(); at_it != bb_atom_list.end(); at_it++){
				atom_shared_ptr mot_at = this_motif->get_bb_atom(*at_it, iter->first);
				if (mot_at->isSet() && mot_at->isActive() && mot_at->get_occupancy() > 0){
					atom_shared_ptr pose_at = pose_->get_bb_atom(*at_it, iter->second);
					if (pose_at->isActive() && pose_at->isSet()){
						atom_mapping[mot_at] = pose_at;

					}
				}
			}

		}
	}


}




/*
bool incr_int_vector(int_vector& vec, const int_vector& upper_limits){
	const int num_ele = vec.size();
	//const int pose_len = rescount;

	for (int i = 0 ; i < num_ele; i++){
		if (vec[i] < upper_limits[i] -1){
			vec[i]++;
			break;
		}
		else if (i < num_ele-1){
			//reset if not last position
			vec[i] = 0;
		}
		else {
			// is last position and overflow
			vec[i]++;
			return false;
		}

	}

	bool overflow=true;
	for (int i = 0 ; i < num_ele; i++){
		overflow = overflow && (vec[i] >= upper_limits[i]);
	}
	if (overflow) return false;


	return true;
}
*/


bool as_geom_pot::assign_motif_by_first_res(POSE::pose_shared_ptr pose_, const int_vector& vec){
	const int mot_chain_count = this->motif_pdb->get_chain_count();
	int pept_mot_chain_count = mot_chain_count;
    for (int ele = 0; ele < this->motif_pdb->get_chain_count(); ele++){
            chain_shared_ptr currChain = motif_pdb->get_chain(ele);
            if (!currChain->isPeptide()){
            	pept_mot_chain_count--;
            }
    }
	const int ele_count = vec.size();
	if (pept_mot_chain_count != ele_count){
		cerr << "as_geom_pot::assign_motif_by_first_res: ERROR: vector size must match element count in motif" << endl;
		return false;
	}
    bool overall_success = true;

    motif_target_res_map_vector.clear();

    const int num_prot_res = pose_->get_residue_count();
    bool_vector assigned_vec(num_prot_res, false);
    if (loop_mask.size() == assigned_vec.size()){
            assigned_vec = loop_mask;
    }


	int pept_ele = 0;
    for (int ele = 0; ele < mot_chain_count; ele++){
            chain_shared_ptr currChain = motif_pdb->get_chain(ele);
            if (currChain->isPeptide()){
                    const int first_res = currChain->get_first_internal_residue_index();
                    const int last_res = currChain->get_last_internal_residue_index();

                    string fr_str = motif_pdb->get_residue(first_res)->get_pdb_residue_index();
                    string lr_str = motif_pdb->get_residue(last_res)->get_pdb_residue_index();

                    trim(fr_str);
                    trim(lr_str);


                    const int fr_num_mot = lexical_cast<int>(fr_str);

                    const int span = lexical_cast<int>(lr_str)
                                            - fr_num_mot + 1;



                    bool is_assigned = false;

                    const int try_pos = vec[pept_ele];

                    if (try_pos+span > pose_->get_residue_count()){
                    	return false;
                    }

                    bool assignOK = true;

                    int_int_map assign_map;

                    chain_shared_ptr fr_chain = pose_->get_residue(try_pos)->get_chain();

                    for (int i = first_res; i <= last_res; i++){
                    	string this_r_str = motif_pdb->get_residue(i)->get_pdb_residue_index();
                    	trim(this_r_str);
                    	const int this_r = lexical_cast<int>(this_r_str);
                    	const int rel_this_r = this_r - fr_num_mot;
                    	if (assigned_vec[try_pos + rel_this_r] == true
                    			|| fr_chain != pose_->get_residue(try_pos + rel_this_r)->get_chain()){
                    		assignOK = false;
                    		break;
                    	}
                    	else {
                    		assign_map[i] = try_pos + rel_this_r;
                    	}
                    }

                    if (assignOK == true){
                    	is_assigned = true;
                    	motif_target_res_map_vector.push_back(assign_map);
                    	int_int_map::const_iterator iter;
                    	for (iter = assign_map.begin(); iter != assign_map.end(); iter++){
                    		assigned_vec[iter->second] = true;
                    	}
                    }



                    if (is_assigned == false){
                            overall_success = false;
                    }




                    pept_ele++;
            }
    }

    if (overall_success == true && motif_target_res_map_vector.size() != (unsigned int)ele_count){
    	cerr << "as_geom_pot::assign_motif_by_first_res: ERROR: oh bugger something weird has gone wrong" << endl;
    	PRINT_EXPR(motif_target_res_map_vector.size());
    	PRINT_EXPR(ele_count);
    	cout << endl << "hmmm\t";
		for (int en = 0; en < ele_count; en++){
			cout << vec[en] << "\t";
		}
		cout << endl;
    	overall_success = false;
    }

    return overall_success;
}

//returns RMSD
double as_geom_pot::do_exaustive_bb_motif_search(POSE::pose_shared_ptr pose_){
	// TODO: do this function recursively somehow
	this->random_assign_motif(pose_);
	double best_rmsd = std::numeric_limits<double>::max();
	int_int_map_vector best_motif_target_res_map_vector;

	const int num_ele = motif_target_res_map_vector.size();
	const int pose_len = pose_->get_residue_count();

	int_vector res_pos(num_ele,0);
	int_vector best_res_pos(num_ele,0);
	const int_vector upper_limit(num_ele,pose_len);
	const int comparisons = UTILS::int_pow(pose_len, num_ele);
	PRINT_EXPR(comparisons);
	for (int i = 0; i < comparisons; i++){
		/*
		PRINT_EXPR(i);
		cout << "TESTING\t" << i << "\t";
		for (int en = 0; en < num_ele; en++){
			cout << res_pos[en] << "\t";
		}
		cout << endl;
		*/


		if (assign_motif_by_first_res(pose_, res_pos)){
			//assigned ok
			const double this_rmsd = this->get_bb_rmsd(pose_);
			if (this_rmsd < best_rmsd){
				best_rmsd = this_rmsd;
				best_motif_target_res_map_vector = motif_target_res_map_vector;
				best_res_pos = res_pos;
			}
		}

		UTILS::incr_vector(res_pos, upper_limit);

	}
	motif_target_res_map_vector = best_motif_target_res_map_vector;


    if (motif_target_res_map_vector.size() != (unsigned int)num_ele){
    	cerr << "as_geom_pot::do_exaustive_bb_motif_search: ERROR: oh bugger something weird has gone wrong" << endl;
    	PRINT_EXPR(motif_target_res_map_vector.size());
    	PRINT_EXPR(num_ele);
    	cout << endl << "hmmm\t";
		for (int en = 0; en < num_ele; en++){
			cout << best_res_pos[en] << "\t";
		}
		cout << endl;
    }

	return best_rmsd;
}

double as_geom_pot::get_bb_rmsd(POSE::pose_shared_ptr pose_) const{
	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;
	this->get_bb_atom_mapping(pose_, atom_mapping);
	return POSE_UTILS::get_rmsd(atom_mapping);
}

pose_shared_ptr as_geom_pot::get_aligned_motif(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_) const{

	pose_shared_ptr this_motif = motif_pdb->clone();

	POSE::pose_shared_ptr pose_ = pose_meta_->get_pose();

	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;

	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		this->get_ca_atom_mapping(pose_, atom_mapping, this_motif);
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta){
		this->get_bb_atom_mapping(pose_, atom_mapping, this_motif);
	}

	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> inv_atom_mapping;
	for (std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr>::const_iterator it = atom_mapping.begin(); it != atom_mapping.end(); it++){
		inv_atom_mapping[it->second] = it->first;
	}

	POSE_UTILS::get_rmsd_superpose(this_motif, inv_atom_mapping);

	return this_motif;
}


double as_geom_pot::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	//double energy = 0;

	POSE::pose_shared_ptr pose_ = pose_meta_->get_pose();

	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;

	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		this->get_ca_atom_mapping(pose_, atom_mapping);
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta){
		this->get_bb_atom_mapping(pose_, atom_mapping);
	}


	const double total_energy = static_cast<double>(atom_mapping.size()) * POSE_UTILS::get_msd(atom_mapping);

	return energies_map.add_energy_component(name_vector[0], total_energy);
	//return energy;
}

double as_geom_pot::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	POSE::pose_shared_ptr pose_ = pose_meta_->get_pose();

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();

	const double weight = energies_map.get_weight(this->name_vector[0]);

	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;

	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		this->get_ca_atom_mapping(pose_, atom_mapping);
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta){
		this->get_bb_atom_mapping(pose_, atom_mapping);
	}

	const double num_atoms = static_cast<double>(atom_mapping.size());
	const double total_energy = num_atoms * POSE_UTILS::get_msd(atom_mapping);

	for (std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr>::iterator at_it = atom_mapping.begin(); at_it != atom_mapping.end(); at_it++){
		atom_shared_ptr currAtom = at_it->second;
		const vector3d init_vec = currAtom->get_coords();
		vector3d vec = init_vec;
		const int index = currAtom->get_seq_num();
		for (int i = 0; i < 3; i++){

			double temp = init_vec[i] + grad_h;
			dummy_funct(temp);
			const double hh =  temp - vec[i];

			vec[i] = init_vec[i];
			vec[i] += hh;
			currAtom->set_coords(vec);
			const double new_val = num_atoms * POSE_UTILS::get_msd(atom_mapping);

			vec[i] = init_vec[i];
			vec[i] -= hh;
			currAtom->set_coords(vec);
			const double new_low_val = num_atoms * POSE_UTILS::get_msd(atom_mapping);

			currAtom->set_coords(init_vec);
			const double diff = new_val - new_low_val;
			const double this_grad = diff / (2.0 * hh);
			(grad[index])[i] += weight * this_grad;

		}
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
	//return energy;

}


void as_geom_pot::dummy_funct(double val) const{
	//double temp = val;
}


bool as_geom_pot::init(){

	if (PRODART::ENV::is_set("as_geom_pot:motif_pdb") == true){
		const string db_path = PRODART::ENV::get_option_value<string>("as_geom_pot:motif_pdb");
		std::ifstream inputdb(db_path.c_str(), ios::in);
		if (inputdb.is_open()){
			this->load_data(inputdb);
			inputdb.close();
			std::cout << "as_geom_pot: loaded motif_pdb file: " << db_path << endl;
			motif_loaded = true;

		}
		else {
			std::cerr << "as_geom_pot: ERROR: can not open motif_pdb file: " << db_path << endl;
			return false;
		}
	}

	if (PRODART::ENV::is_set("as_geom_pot:mapping") == true){
		const string db_path = PRODART::ENV::get_option_value<string>("as_geom_pot:mapping");
		std::ifstream inputdb(db_path.c_str(), ios::in);
		if (inputdb.is_open()){
			this->load_mapping(inputdb);
			inputdb.close();
			std::cout << "as_geom_pot: loaded mapping file: " << db_path << endl;
			motif_loaded = true;

		}
		else {
			std::cerr << "as_geom_pot: ERROR: can not open mapping file: " << db_path << endl;
			return false;
		}
	}

	return true;
}

std::istream& as_geom_pot::load_data( std::istream& input ){

	motif_pdb->loadPdb(input);

	return input;
}

std::istream& as_geom_pot::load_mapping( std::istream& input ){

	this->using_fixed_mapping = true;
	fixed_mapping.clear();
    string lineStr;

    //long length = 0;//, lineNum = 0 ;
	//parse it here
	while ( !input.eof() ) {
		getline(input, lineStr);

		//length = lineStr.length();

		string_vector SplitVec;
		split( SplitVec, lineStr, is_any_of("\t ") );

		if (SplitVec.size() == 4){
			if ( SplitVec[0].substr(0,1).compare("#") != 0){
				const string ref_resnum = SplitVec[0];
				char ref_chain = lexical_cast<char>(SplitVec[1]);
				if (ref_chain == '-') ref_chain = ' ';
				const string resnum = SplitVec[2];
				char tm_chain = lexical_cast<char>(SplitVec[3]);
				if (tm_chain == '-') tm_chain = ' ';

				string_char_string_char_tuple tup(ref_resnum, ref_chain, resnum, tm_chain);
				fixed_mapping.push_back(tup);

				cout << "as_geom_pot: load_mapping: fixed mapping: " << ref_resnum << "\t"
						<< ref_chain << "\t"
						<< resnum << "\t"
						<< tm_chain << "\t"
						<< "\n";


			}

		}

	}


	return input;
}

bool as_geom_pot::apply_fixed_mapping(POSE::pose_shared_ptr pose_){
	if (!using_fixed_mapping) return false;
	motif_target_res_map_vector.clear();


	int_int_map assign_map;
	for (unsigned int i = 0; i < fixed_mapping.size(); i++){
		string mot_res = fixed_mapping[i].get<0>();
		char mot_ch = fixed_mapping[i].get<1>();
		string pose_res = fixed_mapping[i].get<2>();
		char pose_ch = fixed_mapping[i].get<3>();

		residue_shared_ptr mot_ptr = motif_pdb->get_residue(mot_res, mot_ch);
		residue_shared_ptr pose_ptr = pose_->get_residue(pose_res, pose_ch);

		if (mot_ptr && pose_ptr){
			assign_map[mot_ptr->get_internal_residue_index()] = pose_ptr->get_internal_residue_index();
		}
		else {
			cerr << "as_geom_pot: apply_fixed_mapping: ERROR: could not find mapped atoms" << endl;
			return false;
		}


	}
	motif_target_res_map_vector.push_back(assign_map);


	return true;
}




bool as_geom_pot::random_assign_motif(POSE::pose_shared_ptr pose_){

	//cout << "random_assign_motif" << endl;

	if (using_fixed_mapping){
		return apply_fixed_mapping(pose_);
	}

	bool overall_success = false;

	const int max_overall_tries = 100;
	int overall_tries = 0;

	while (overall_success == false && overall_tries < max_overall_tries){

		overall_success = true;

		motif_target_res_map_vector.clear();

		const int num_prot_res = pose_->get_residue_count();
		bool_vector assigned_vec(num_prot_res, false);
		if (loop_mask.size() == assigned_vec.size()){
			assigned_vec = loop_mask;
		}

		const int mot_chain_count = this->motif_pdb->get_chain_count();

		for (int ele = 0; ele < mot_chain_count; ele++){
			chain_shared_ptr currChain = motif_pdb->get_chain(ele);
			if (currChain->isPeptide()){
				const int first_res = currChain->get_first_internal_residue_index();
				const int last_res = currChain->get_last_internal_residue_index();

				string fr_str = motif_pdb->get_residue(first_res)->get_pdb_residue_index();
				string lr_str = motif_pdb->get_residue(last_res)->get_pdb_residue_index();

				trim(fr_str);
				trim(lr_str);

				const int fr_num_mot = lexical_cast<int>(fr_str);

				const int span = lexical_cast<int>(lr_str)
							- fr_num_mot + 1;

				bool is_assigned = false;

				const int max_ele_tries = 500;
				int ele_tries = 0;

				// assign element
				while (is_assigned == false && ele_tries < max_ele_tries){

					const int try_pos = rand_num->randInt(num_prot_res - span);
					bool assignOK = true;

					int_int_map assign_map;

					chain_shared_ptr fr_chain = pose_->get_residue(try_pos)->get_chain();

					for (int i = first_res; i <= last_res; i++){
						string this_r_str = motif_pdb->get_residue(i)->get_pdb_residue_index();
						trim(this_r_str);
						const int this_r = lexical_cast<int>(this_r_str);
						const int rel_this_r = this_r - fr_num_mot;
						if (assigned_vec[try_pos + rel_this_r] == true
								|| fr_chain != pose_->get_residue(try_pos + rel_this_r)->get_chain()){
							assignOK = false;
							break;
						}
						else {
							assign_map[i] = try_pos + rel_this_r;
						}
					}

					if (assignOK == true){
						is_assigned = true;
						motif_target_res_map_vector.push_back(assign_map);
						int_int_map::const_iterator iter;
						for (iter = assign_map.begin(); iter != assign_map.end(); iter++){
							assigned_vec[iter->second] = true;
						}
					}
					ele_tries++;
				}

				if (is_assigned == false){
					overall_success = false;
				}

			}
		}
		overall_tries++;
	}

	return overall_success;
}


void as_geom_pot::print_assign_info(std::ostream& output) const{

	for (unsigned int ele = 0; ele < motif_target_res_map_vector.size(); ele++){
		int_int_map::const_iterator iter;
		for (iter = motif_target_res_map_vector[ele].begin(); iter != motif_target_res_map_vector[ele].end(); iter++ ){

			output << ele << "\t"
					<< iter->first << "\t"
					<< iter->second << "\t"
					<< endl;


		}

	}

}


void as_geom_pot::temp_store_mapping(){
	temp_bk_motif_target_res_map_vector = motif_target_res_map_vector; //atom_mapping;
}

void as_geom_pot::temp_restore_mapping(){
	motif_target_res_map_vector = temp_bk_motif_target_res_map_vector;
}

void as_geom_pot::move_store_mapping(){
	move_bk_motif_target_res_map_vector = motif_target_res_map_vector;
}
void as_geom_pot::move_restore_mapping(){
	motif_target_res_map_vector = move_bk_motif_target_res_map_vector;
}

bool as_geom_pot::move_element(POSE::pose_shared_ptr pose_,
		const int ele,
		const int max_move){



	const int num_ele = motif_target_res_map_vector.size();

	if (ele >= num_ele){
		return false;
	}

	this->temp_store_mapping();
	const int num_prot_res = pose_->get_residue_count();

	bool_vector assigned_vec(num_prot_res, false);
	if (loop_mask.size() == assigned_vec.size()){
		assigned_vec = loop_mask;
	}
	for (int i = 0; i < num_ele; i++ ){
		if (i != ele){
			for (int_int_map::const_iterator at_it = motif_target_res_map_vector[i].begin(); at_it != motif_target_res_map_vector[i].end(); at_it++){
				assigned_vec[at_it->second] = true;
			}
		}
	}

	bool is_assigned = false;
	const int max_ele_tries = 10;
	int ele_tries = 0;

	// assign element
	while (is_assigned == false && ele_tries < max_ele_tries){


		int move = 1 + rand_num->randInt(max_move - 1);

		if (rand_num->randInt(1) == 0){
			move = -move;
		}

		int_int_map new_assignment;

		bool assignOK = true;
		for (int_int_map::const_iterator at_it = motif_target_res_map_vector[ele].begin(); at_it != motif_target_res_map_vector[ele].end(); at_it++){
			const int new_pos = at_it->second + move;
			if (new_pos < 0 || new_pos >= num_prot_res){
				assignOK = false;
				break;
			}
			else if (assigned_vec[new_pos] == true ){
				assignOK = false;
				break;
			}
			else {
				//at_it->second += move;
				new_assignment[at_it->first] = at_it->second + move;
			}
		}

		chain_shared_ptr currChain = pose_->get_residue(new_assignment.begin()->second)->get_chain();
		for (int_int_map::const_iterator at_it = new_assignment.begin(); at_it != new_assignment.end(); at_it++){
			if (pose_->get_residue(at_it->second)->get_chain() != currChain){
				assignOK = false;
				break;
			}
		}

		if (assignOK == true){
			is_assigned = true;
			motif_target_res_map_vector[ele] = new_assignment;
			//cout << "move_element.  ele: " << ele << " max_move: " << max_move << " move: " << move << " tries: " << ele_tries << endl;
		}

		ele_tries++;
	}

	return is_assigned;
}


bool as_geom_pot::make_move(POSE::pose_shared_ptr pose_){


	if (using_fixed_mapping){
		return apply_fixed_mapping(pose_);
	}

	const double choice = rand_num->rand();
	const int max_tries = 100;
	int tries = 0;
	bool result = false;

	while (result == false && tries <= max_tries){
		if (choice < 0.5){
			//small move
			const int num_ele = motif_target_res_map_vector.size();
			const int ele = rand_num->randInt(num_ele-1);
			result = this->move_element(pose_, ele, 1);
		}
		else if (choice < 0.95){
			//large move
			const int num_ele = motif_target_res_map_vector.size();
			const int ele = rand_num->randInt(num_ele-1);
			result = this->move_element(pose_, ele, 5);
		}
		else if (choice < 0.975){
			//large move
			const int num_ele = motif_target_res_map_vector.size();
			const int ele = rand_num->randInt(num_ele-1);
			result = this->move_element(pose_, ele, pose_->get_residue_count() / 4);
		}
		else {
			result = this->random_assign_motif(pose_);
		}
		tries++;
	}

	return result;
}

void as_geom_pot::paint_motif_occupancy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_, const double base_value ) {

	pose_shared_ptr pose_ = pose_meta_->get_pose();

	//const int resCount = pose_->get_residue_count();
	const int prot_size = pose_->get_all_atom_count();

	for (int i = 0; i < prot_size; i++){
		atom_shared_ptr currAtom = pose_->get_atom(i);

		if (currAtom){
			currAtom->set_occupancy( base_value );
		}

	}

	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;

	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		this->get_ca_atom_mapping(pose_, atom_mapping);
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta){
		this->get_bb_atom_mapping(pose_, atom_mapping);
	}

	for (std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr>::const_iterator it = atom_mapping.begin(); it != atom_mapping.end(); it++){
		it->second->set_occupancy(base_value + 1.0);
	}

	/*
	const int motif_size = motif_coors_vec.size();

	for (int i = 0; i< motif_size; i++){
		Atom* currAtom = protein.getAtomByAtomIndex(atom_mapping[i]);
		const int ele_index = motif_coors_vec[i].element_index;
		if (currAtom != 0){
			currAtom->set_occupancy( base_value + 1.0 + ele_index );
		}
	}
	*/
}


void as_geom_pot::transfer_sidechains(const PRODART::POSE::META::bb_pose_meta_shared_ptr pose_meta_) const{
	pose_shared_ptr pose_ = pose_meta_->get_pose();
	pose_shared_ptr this_motif = this->get_aligned_motif(pose_meta_);

	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;
	this->get_bb_atom_mapping(pose_, atom_mapping, this_motif);
	std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr>::iterator it;

	for (it = atom_mapping.begin(); it != atom_mapping.end(); it++){
		if (it->first->get_type()== atom_type("CB") && it->second->get_type()== atom_type("CB")){
			sidechain_shared_ptr pose_sc = it->second->get_residue()->get_sidechain();
			sidechain_shared_ptr mot_sc = it->first->get_residue()->get_sidechain();
			cout << "copying: " << it->first->get_residue()->get_pdb_residue_index() << "\t to " << it->second->get_residue()->get_pdb_residue_index()
					<< endl;

			pose_sc->copy_sidechain(mot_sc);
			it->second->get_residue()->set_type(it->first->get_residue()->get_type());

			stringstream sstrm (stringstream::in | stringstream::out);

			sstrm << "as_geom_pot: motif_to_scaffold_mapping:\t"
					<<  it->first->get_residue()->get_pdb_residue_index()
					<< "\t" << it->first->get_residue()->get_chain()->getChainID()
					<< "\t" << it->first->get_residue()->get_type().get_pdb_formatted_label()
					<< "\t->\t"
					<< "\t" << it->second->get_residue()->get_pdb_residue_index()
					<< "\t" << it->second->get_residue()->get_chain()->getChainID()
					<< "\t" << it->second->get_residue()->get_type().get_pdb_formatted_label();
			pose_->add_remark_with_prodart_label(sstrm.str());
		}
	}

	pose_->index();
}

void as_geom_pot::transfer_non_peptides(const PRODART::POSE::META::bb_pose_meta_shared_ptr pose_meta_) const{
	pose_shared_ptr pose_ = pose_meta_->get_pose();
	pose_shared_ptr this_motif = this->get_aligned_motif(pose_meta_);

	const int num_chains = this_motif->get_chain_count();
	for (int i = 0; i < num_chains; i++){
		chain_shared_ptr currChain = this_motif->get_chain(i);
		if (!currChain->isPeptide()){
			cout << "copying chain " << currChain->getChainID() << endl;
			pose_->add_duplicated_chain(currChain);
		}
	}

	pose_->index();
}


}
}
}
