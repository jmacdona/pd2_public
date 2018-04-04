/*
 * ca2main_protocols.cpp
 *
 *  Created on: 24 Sep 2010
 *      Author: jmacdona
 */
#include "ca2main_protocols.h"
#include "protocols.h"

using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::SIM;

namespace PRODART {
namespace PROTOCOLS{
namespace CA2MAIN{

vector3d new_ca_pos(vector3d c1,
		vector3d c2,
		vector3d c3,
		vector3d c4){

	const double dih = dihedral(c1, c2, c3, c4);
	const double ang = angle(c2, c3, c4);
	const double bond_len = 3.8;
	vector3d new_ca = dihedralEnd(c2, c3, c4, bond_len, ang, dih);
	return new_ca;
}

bool add_two_end_ca_atoms_residues(PRODART::POSE::pose_shared_ptr protein){

	const int num_chains = protein->get_chain_count();

	for (int i = 0; i < num_chains; i++  ){
		chain_shared_ptr currChain = protein->get_chain(i);
		const int orig_len = currChain->length();
		//cout << "a " << currChain->get_first_internal_residue_index() << "\t" << currChain->get_last_internal_residue_index() << endl;
		if (orig_len > 0 && currChain->isPeptide()){
			//cout << "adding to chain: " << currChain->getChainID() << endl;
			residue_shared_ptr  n2 = protein->add_nterm_residue(residue_type("GLY"), i);
			//cout << "a " << currChain->get_first_internal_residue_index() << "\t" << currChain->get_last_internal_residue_index() << endl;
			residue_shared_ptr  n1 = protein->add_nterm_residue(residue_type("GLY"), i);
			//cout << "b " << currChain->get_first_internal_residue_index() << "\t" << currChain->get_last_internal_residue_index() << endl;
			residue_shared_ptr  c1 = protein->add_cterm_residue(residue_type("GLY") , i);
			//cout << "c " << currChain->get_first_internal_residue_index() << "\t" << currChain->get_last_internal_residue_index() << endl;
			residue_shared_ptr  c2 = protein->add_cterm_residue(residue_type("GLY") , i);
			//cout << "d " << currChain->get_first_internal_residue_index() << "\t" << currChain->get_last_internal_residue_index() << endl;

			if (orig_len >= 4 ){
				//cout << "adding n2" << endl;
				for (int i = 2; i <= 5; i++){
					atom_shared_ptr ca_test = currChain->get_ca(i);
					if (!ca_test){
						cout << "ERROR: problem with atom " << i << endl;
					}

					//cout << currChain->get_first_internal_residue_index() << "\t" << currChain->get_last_internal_residue_index() << endl;

				}
				const vector3d new_ca_coords_n2 = new_ca_pos(currChain->get_ca_pos(5),
						currChain->get_ca_pos(4),
						currChain->get_ca_pos(3),
						currChain->get_ca_pos(2));
				protein->add_new_atom(new_ca_coords_n2, atom_type("CA"), n2->get_internal_residue_index());
				//cout << "adding n1" << endl;
				const vector3d new_ca_coords_n1 = new_ca_pos(currChain->get_ca_pos(4),
						currChain->get_ca_pos(3),
						currChain->get_ca_pos(2),
						currChain->get_ca_pos(1));
				protein->add_new_atom(new_ca_coords_n1, atom_type("CA"), n1->get_internal_residue_index());

				//cout << "adding c1" << endl;
				const int len = currChain->length();
				const vector3d new_ca_coords_c1 = new_ca_pos(currChain->get_ca_pos(len - 6),
						currChain->get_ca_pos(len - 5),
						currChain->get_ca_pos(len - 4),
						currChain->get_ca_pos(len - 3));
				protein->add_new_atom(new_ca_coords_c1, atom_type("CA"), c1->get_internal_residue_index());
				//cout << "adding c2" << endl;
				const vector3d new_ca_coords_c2 = new_ca_pos(currChain->get_ca_pos(len - 5),
						currChain->get_ca_pos(len - 4),
						currChain->get_ca_pos(len - 3),
						currChain->get_ca_pos(len - 2));
				protein->add_new_atom(new_ca_coords_c2, atom_type("CA"), c2->get_internal_residue_index());
				protein->index();

			}
			else {
				// TODO some other add method for short chains

				protein->index();
			}

			// TODO add step to relax newly added terminal residues


		}
	}

	return true;
}

bool delete_two_end_residues(PRODART::POSE::pose_shared_ptr protein){
	const int num_chains = protein->get_chain_count();

	for (int i = 0; i < num_chains; i++  ){
		chain_shared_ptr currChain = protein->get_chain(i);
		if (currChain->length() > 4 && currChain->isPeptide()){
			protein->delete_residue(currChain->get_first_internal_residue_index());
			protein->delete_residue(currChain->get_first_internal_residue_index());
			protein->delete_residue(currChain->get_last_internal_residue_index());
			protein->delete_residue(currChain->get_last_internal_residue_index());
		}

	}

	protein->index();
	return true;
}

bool simple_ca2main(PRODART::POSE::pose_shared_ptr protein){

	cout << "PRODART::PROTOCOLS::CA2MAIN::simple_ca2main : running..." << endl;

	PRODART::POSE::BB_BUILDER::const_backbone_builder_shared_ptr bb_builder = PRODART::POSE::BB_BUILDER::backbone_builder::Instance();
	cout << "simple_ca2main: pose size: " << protein->get_residue_count() << endl;
	cout << "simple_ca2main: number of chains: " << protein->get_chain_count() << endl;
	add_two_end_ca_atoms_residues(protein);
	//cout << "pose size " << protein->get_residue_count() << endl;
	PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(protein);
	bb_builder->buildBackbone(meta_data);
	delete_two_end_residues(protein);
	meta_data->inactivate_pseudo_NO();

	cout << "PRODART::PROTOCOLS::CA2MAIN::simple_ca2main : done." << endl;

	return true;
}

bool simple_ca2main(PRODART::POSE::pose_shared_ptr protein, const bool_vector& loop_mask){

	PRODART::POSE::META::ca_pose_meta_shared_ptr ca_meta_data = PRODART::POSE::META::new_ca_pose_meta(protein);
	PRODART::POSE::BB_BUILDER::const_backbone_builder_shared_ptr bb_builder = PRODART::POSE::BB_BUILDER::backbone_builder::Instance();

	frag4_vector& fragments = ca_meta_data->get_fragments();

	frag4_vector::iterator iter;
	for (iter = fragments.begin(); iter != fragments.end(); iter++){
		if (loop_mask[iter->CA_1_residue_number] || loop_mask[iter->CA_1_residue_number+1]){
			bb_builder->buildBackbone_single(ca_meta_data, *iter, loop_mask[iter->CA_1_residue_number]);
		}
	}
	ca_meta_data->inactivate_pseudo_NO();
	MISC::inactivate_unset_GLY_CBs(protein);
	return true;
}

// new ca2main with no CA refinement, minimisation or optimisation
bool new_ca2main(PRODART::POSE::pose_shared_ptr protein){
	cout << "PRODART::PROTOCOLS::CA2MAIN::new_ca2main : running..." << endl;

	cout << "new_ca2main: pose size: " << protein->get_residue_count() << endl;
	cout << "new_ca2main: number of chains: " << protein->get_chain_count() << endl;
	add_two_end_ca_atoms_residues(protein);
	add_two_end_ca_atoms_residues(protein);
	BB_BUILDER::alphabet_bb_builder_manager::Instance()->build_backbone(protein);
	delete_two_end_residues(protein);
	delete_two_end_residues(protein);
	PRODART::POSE::META::ca_pose_meta::inactivate_pseudo_NO(protein);
	cout << "PRODART::PROTOCOLS::CA2MAIN::new_ca2main : done." << endl;

	return true;
}

bool new_ca2main(PRODART::POSE::pose_shared_ptr protein, const bool_vector& loop_mask){

	cerr << "\nnew_ca2main: warning: this routine MAY BE BUGGY but seems to be ok so far\n" << endl;

	int_vector starts, ends;
	bool inEle = false;
	bool last_val = false;
	cout << "new_ca2main: loop_mask: ";
	for (int i = 0 ; i < loop_mask.size(); i++){
		cout << loop_mask[i];
	}
	cout << endl;
	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if ((last_val != loop_mask[i]) && loop_mask[i]){
			//start of element
			inEle = true;
			// push back one before moved CA to cover previous peptide bond
			starts.push_back(i > 0 ? i-1 : 0);
			cerr << "start " << starts.back() << endl;
		}
		else if ((last_val != loop_mask[i]) && !loop_mask[i]){
			// gone past end of element
			inEle = false;
			ends.push_back(i < loop_mask.size() ? i : loop_mask.size());
			cerr << "end " << ends.back() << endl;
		}

		last_val = loop_mask[i];
	}
	if (inEle){
		ends.push_back(loop_mask.size()-1);
		cerr << "end " << ends.back() << endl;
	}

	if (starts.size() != ends.size()){
		cerr << "new_ca2main: ERROR: some kind of error parsing loop_mask" << endl;
		PRINT_EXPR(starts.size());
		PRINT_EXPR(ends.size());
		//PRINT_EXPR(loop_mask);
		return false;
	}

	PRINT_EXPR(loop_mask.size());
	PRINT_EXPR(protein->get_residue_count());

	for (unsigned int i = 0; i < starts.size(); i++){
		cerr << "internal: " << starts[i] << "\t"
				<< ends[i] << "\t"
				<< endl;
		cerr << "pdb: " << protein->get_residue(starts[i])->get_trimmed_pdb_residue_index() << "\t"
				<< protein->get_residue(ends[i])->get_trimmed_pdb_residue_index() << "\t"
				<< endl;
		BB_BUILDER::alphabet_bb_builder_manager::Instance()->build_peptide_bonds(starts[i], ends[i], protein);

	}

	//cerr << "\nnew_ca2main: warning: this routine MAY BE BUGGY but seems to be ok so far\n" << endl;

	PRODART::POSE::META::ca_pose_meta::inactivate_pseudo_NO(protein);
	return true;
}

// simple ca2main with no CA refinement, minimisation or optimisation
bool simple_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps,
		const std::string bb_min_potentials_weight_set){
	simple_ca2main(protein);
	bb_minimise(protein, bb_min_potentials_weight_set, min_steps);
	return true;
}


// simple ca2main with minimisation
bool simple_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot){
	simple_ca2main(protein, loop_mask);
	bb_minimise(protein, loop_mask, bb_min_pot, min_steps);
	return true;
}

//  ca2main with minimisation
bool new_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot){
	new_ca2main(protein, loop_mask);
	bb_minimise(protein, loop_mask, bb_min_pot, min_steps);
	return true;
}

// alphabet ca2main with minimisation
bool new_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps,
		const std::string bb_min_pot){
	new_ca2main(protein);
	bb_minimise(protein, bb_min_pot, min_steps);
	return true;
}

bool new_ca2main_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps,
		const std::string bb_min_pot){
	new_ca2main(protein);
	bb_minimise_fixed_ca(protein, bb_min_pot, min_steps);
	return true;
}

bool new_ca2main_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		const std::string bb_min_pot){
	new_ca2main(protein, loop_mask);
	bb_minimise_fixed_ca(protein, loop_mask, bb_min_pot, min_steps);
	return true;
}

bool simple_ca2main_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		const std::string bb_min_pot){
	simple_ca2main(protein, loop_mask);
	bb_minimise_fixed_ca(protein, loop_mask, bb_min_pot, min_steps);
	return true;
}

bool new_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		const std::string bb_min_potentials_weight_set){
	new_ca2main(protein, loop_mask);
	bb_minimise(protein, loop_mask, bb_min_potentials_weight_set, min_steps);
	return true;
}

// simple ca2main with no CA refinement, minimisation or optimisation
bool simple_ca2main_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps,
		const std::string bb_min_potentials_weight_set){
	simple_ca2main(protein, loop_mask);
	bb_minimise(protein, loop_mask, bb_min_potentials_weight_set, min_steps);
	return true;
}

bool bb_minimise(PRODART::POSE::pose_shared_ptr protein,
		const std::string pot_preset,
		const unsigned long min_steps){
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container(pot_preset);
	potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);

	minimiser_interface_shared_ptr min_sim =  new_minimiser(bb_min_pot);
	min_sim->set_steps(min_steps);

	const bool result = min_sim->make_move(bb_meta);

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);


	return result;
}

PRODART::POSE::atom_shared_ptr_vector get_loop_moved_atoms(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask){
	atom_shared_ptr_vector moved_atoms;


	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){

			if (i > 0){
				if (protein->get_bb_atom(C, i-1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(C, i-1));
				if (protein->get_bb_atom(O, i-1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(O, i-1));
			}

			if (protein->get_bb_atom(N, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(N, i));
			if (protein->get_bb_atom(POSE::CA, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(POSE::CA, i));
			if (protein->get_bb_atom(C, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(C, i));
			if (protein->get_bb_atom(CB, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(CB, i));
			if (protein->get_bb_atom(O, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(O, i));
			if (protein->get_bb_atom(H, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(H, i));
			if (i+1 < loop_mask.size()){
				if (protein->get_bb_atom(N, i+1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(N, i+1));
				if (protein->get_bb_atom(H, i+1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(H, i+1));
			}
		}
	}

	return moved_atoms;

}

PRODART::POSE::atom_shared_ptr_vector get_loop_moved_N_CA_C_O_atoms(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask){
	atom_shared_ptr_vector moved_atoms;

	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){

			if (i > 0){
				if (protein->get_bb_atom(C, i-1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(C, i-1));
				if (protein->get_bb_atom(O, i-1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(O, i-1));
			}

			if (protein->get_bb_atom(N, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(N, i));
			if (protein->get_bb_atom(POSE::CA, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(POSE::CA, i));
			if (protein->get_bb_atom(C, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(C, i));
			//if (protein->get_bb_atom(CB, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(CB, i));
			if (protein->get_bb_atom(O, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(O, i));
			//if (protein->get_bb_atom(H, i)->isActive()) moved_atoms.push_back(protein->get_bb_atom(H, i));
			if (i+1 < loop_mask.size()){
				if (protein->get_bb_atom(N, i+1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(N, i+1));
				//if (protein->get_bb_atom(H, i+1)->isActive()) moved_atoms.push_back(protein->get_bb_atom(H, i+1));
			}
		}
	}

	return moved_atoms;
}

bool_vector make_residue_loop_mask(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum){
	bool_vector loop_mask(protein->get_residue_count(), false);
	if (end_resnum >= loop_mask.size()){
		cerr << "make_residue_loop_mask: ERROR: end_resnum is more than size of array!" << endl;
	}
	residue_shared_ptr_vector::const_iterator iter;
	for (int i = start_resnum; i <= end_resnum; i++){
		loop_mask[i] = true;
	}
	return loop_mask;
}

bool_vector make_residue_loop_mask(const PRODART::POSE::three_state_sec_struct_vector& secs){
	bool_vector loop_mask(secs.size(), false);
	residue_shared_ptr_vector::const_iterator iter;
	for (unsigned int i = 0; i < secs.size(); i++){
		if (secs[i] == ss3_OTHER) loop_mask[i] = true;
	}
	return loop_mask;
}

PRODART::POSE::atom_shared_ptr_vector get_loop_moved_atoms(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum){
	bool_vector loop_mask = make_residue_loop_mask(protein, start_resnum, end_resnum);
	return get_loop_moved_atoms(protein, loop_mask);
}

PRODART::POSE::atom_shared_ptr_vector get_loop_moved_N_CA_C_O_atoms(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum){
	bool_vector loop_mask = make_residue_loop_mask(protein, start_resnum, end_resnum);
	return get_loop_moved_N_CA_C_O_atoms(protein, loop_mask);
}

bool bb_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot,
		const unsigned long min_steps){
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
#ifndef NDEBUG
	cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
	cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
#endif
	potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();
	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	//cout << endl;

	bool_vector atom_selection(protein->get_all_atom_count(), false);
	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){

			if (i > 0){
				if (protein->get_bb_atom(C, i-1)->isActive()) atom_selection[protein->get_bb_atom(C, i-1)->get_seq_num()] = true;
				if (protein->get_bb_atom(O, i-1)->isActive()) atom_selection[protein->get_bb_atom(O, i-1)->get_seq_num()] = true;
			}

			if (protein->get_bb_atom(N, i)->isActive()) atom_selection[protein->get_bb_atom(N, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(POSE::CA, i)->isActive()) atom_selection[protein->get_bb_atom(POSE::CA, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(C, i)->isActive()) atom_selection[protein->get_bb_atom(C, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(CB, i)->isActive()) atom_selection[protein->get_bb_atom(CB, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(O, i)->isActive()) atom_selection[protein->get_bb_atom(O, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(H, i)->isActive()) atom_selection[protein->get_bb_atom(H, i)->get_seq_num()] = true;
			if (i+1 < loop_mask.size()){
				if (protein->get_bb_atom(N, i+1)->isActive()) atom_selection[protein->get_bb_atom(N, i+1)->get_seq_num()] = true;
				if (protein->get_bb_atom(H, i+1)->isActive()) atom_selection[protein->get_bb_atom(H, i+1)->get_seq_num()] = true;
			}
		}
	}
	minimiser_interface_shared_ptr min_sim =  new_minimiser(bb_min_pot, atom_selection);
	min_sim->set_steps(min_steps);

	const bool result = min_sim->make_move(bb_meta);

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	//enrg_map.print_weights(cout);
	enrg_map.print_weighted_components(cout);
	return result;
}

bool bb_minimise(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const std::string pot_preset,
		const unsigned long min_steps){
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
#ifndef NDEBUG
	cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
	cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
#endif
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container(pot_preset);
	potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();
	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	//cout << endl;

	bool_vector atom_selection(protein->get_all_atom_count(), false);
	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){

			if (i > 0){
				if (protein->get_bb_atom(C, i-1)->isActive()) atom_selection[protein->get_bb_atom(C, i-1)->get_seq_num()] = true;
				if (protein->get_bb_atom(O, i-1)->isActive()) atom_selection[protein->get_bb_atom(O, i-1)->get_seq_num()] = true;
			}

			if (protein->get_bb_atom(N, i)->isActive()) atom_selection[protein->get_bb_atom(N, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(POSE::CA, i)->isActive()) atom_selection[protein->get_bb_atom(POSE::CA, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(C, i)->isActive()) atom_selection[protein->get_bb_atom(C, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(CB, i)->isActive()) atom_selection[protein->get_bb_atom(CB, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(O, i)->isActive()) atom_selection[protein->get_bb_atom(O, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(H, i)->isActive()) atom_selection[protein->get_bb_atom(H, i)->get_seq_num()] = true;
			if (i+1 < loop_mask.size()){
				if (protein->get_bb_atom(N, i+1)->isActive()) atom_selection[protein->get_bb_atom(N, i+1)->get_seq_num()] = true;
				if (protein->get_bb_atom(H, i+1)->isActive()) atom_selection[protein->get_bb_atom(H, i+1)->get_seq_num()] = true;
			}
		}
	}
	minimiser_interface_shared_ptr min_sim =  new_minimiser(bb_min_pot, atom_selection);
	min_sim->set_steps(min_steps);

	const bool result = min_sim->make_move(bb_meta);

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	return result;
}


bool bb_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const std::string pot_preset,
		const unsigned long min_steps){
	bool_vector loop_mask = make_residue_loop_mask( protein,
			0,
			protein->get_residue_count()-1);
	return bb_minimise_fixed_ca(protein, loop_mask, pot_preset, min_steps);

}
bool bb_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot,
		const unsigned long min_steps ){
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
#ifndef NDEBUG
	cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
	cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
#endif
	potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();
	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	//cout << endl;

	bool_vector atom_selection(protein->get_all_atom_count(), false);
	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){

			if (i > 0){
				if (protein->get_bb_atom(C, i-1)->isActive()) atom_selection[protein->get_bb_atom(C, i-1)->get_seq_num()] = true;
				if (protein->get_bb_atom(O, i-1)->isActive()) atom_selection[protein->get_bb_atom(O, i-1)->get_seq_num()] = true;
			}

			if (protein->get_bb_atom(N, i)->isActive()) atom_selection[protein->get_bb_atom(N, i)->get_seq_num()] = true;
			//if (protein->get_bb_atom(POSE::CA, i)->isActive()) atom_selection[protein->get_bb_atom(POSE::CA, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(C, i)->isActive()) atom_selection[protein->get_bb_atom(C, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(CB, i)->isActive()) atom_selection[protein->get_bb_atom(CB, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(O, i)->isActive()) atom_selection[protein->get_bb_atom(O, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(H, i)->isActive()) atom_selection[protein->get_bb_atom(H, i)->get_seq_num()] = true;
			if (i+1 < loop_mask.size()){
				if (protein->get_bb_atom(N, i+1)->isActive()) atom_selection[protein->get_bb_atom(N, i+1)->get_seq_num()] = true;
				if (protein->get_bb_atom(H, i+1)->isActive()) atom_selection[protein->get_bb_atom(H, i+1)->get_seq_num()] = true;
			}
		}
	}
	minimiser_interface_shared_ptr min_sim =  new_minimiser(bb_min_pot, atom_selection);
	min_sim->set_steps(min_steps);

	const bool result = min_sim->make_move(bb_meta);

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	//enrg_map.print_weights(cout);
	enrg_map.print_weighted_components(cout);
	return result;
}
bool bb_minimise_fixed_ca(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const std::string pot_preset,
		const unsigned long min_steps){
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container(pot_preset);
	return bb_minimise_fixed_ca(protein, loop_mask, bb_min_pot, min_steps);
}

bool bb_minimise_geom(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long min_steps){

	return bb_minimise(protein, "bb_min_default", min_steps);
	/*
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("bb_min_default");
	potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);

	minimiser_shared_ptr min_sim =  new_minimiser(bb_min_pot);
	min_sim->set_steps(min_steps);

	const bool result = min_sim->make_move(bb_meta);

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);


	return result;
	*/
}

bool bb_minimise_geom(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long min_steps){

	return bb_minimise(protein, loop_mask, "bb_min_default", min_steps);

	/*
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
	cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
	cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("bb_min_default");
	potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();
	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	//cout << endl;

	bool_vector atom_selection(protein->get_all_atom_count(), false);
	for (unsigned int i = 0; i < loop_mask.size(); i++){
		if (loop_mask[i]){

			if (i > 0){
				if (protein->get_bb_atom(C, i-1)->isActive()) atom_selection[protein->get_bb_atom(C, i-1)->get_seq_num()] = true;
				if (protein->get_bb_atom(O, i-1)->isActive()) atom_selection[protein->get_bb_atom(O, i-1)->get_seq_num()] = true;
			}

			if (protein->get_bb_atom(N, i)->isActive()) atom_selection[protein->get_bb_atom(N, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(POSE::CA, i)->isActive()) atom_selection[protein->get_bb_atom(POSE::CA, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(C, i)->isActive()) atom_selection[protein->get_bb_atom(C, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(CB, i)->isActive()) atom_selection[protein->get_bb_atom(CB, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(O, i)->isActive()) atom_selection[protein->get_bb_atom(O, i)->get_seq_num()] = true;
			if (protein->get_bb_atom(H, i)->isActive()) atom_selection[protein->get_bb_atom(H, i)->get_seq_num()] = true;
			if (i+1 < loop_mask.size()){
				if (protein->get_bb_atom(N, i+1)->isActive()) atom_selection[protein->get_bb_atom(N, i+1)->get_seq_num()] = true;
				if (protein->get_bb_atom(H, i+1)->isActive()) atom_selection[protein->get_bb_atom(H, i+1)->get_seq_num()] = true;
			}
		}
	}
	minimiser_shared_ptr min_sim =  new_minimiser(bb_min_pot, atom_selection);
	min_sim->set_steps(min_steps);

	const bool result = min_sim->make_move(bb_meta);

	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	return result;
	*/
}


int_vector_vector get_strands_vector(const PRODART::POSE::three_state_sec_struct_vector& secs){
	if (secs.size() == 0){
		return int_vector_vector();
	}

	int_vector_vector strands;

	//int_vector::const_iterator iter;
	strands.clear();

	three_state_sec_struct last_val = ss3_UNDEF;
	int num_strands = 0;
	int strand_length = 0;


	const int resNum = secs.size();

	for (int i = 0; i < resNum; i++){
		const three_state_sec_struct this_val = secs[i];

		if (last_val != ss3_STRAND && this_val == ss3_STRAND){
			num_strands++;
			strand_length = 1;
			int_vector temp_vec(1,i);
			strands.push_back(temp_vec);

		}
		else if (this_val == ss3_STRAND){
			strands[num_strands-1].push_back(i);
			strand_length++;
		}


		last_val = this_val;
	}

	return strands;
}
int_vector_vector get_helices_vector(const PRODART::POSE::three_state_sec_struct_vector& secs){
	if (secs.size() == 0){
		return int_vector_vector();
	}

	int_vector_vector helices;

	//int_vector::const_iterator iter;
	helices.clear();

	three_state_sec_struct last_val = ss3_UNDEF;
	int num_helices = 0;
	int helix_length = 0;


	const int resNum = secs.size();

	for (int i = 0; i < resNum; i++){
		const three_state_sec_struct this_val = secs[i];

		if (last_val != ss3_HELIX && this_val == ss3_HELIX){
			num_helices++;
			helix_length = 1;
			int_vector temp_vec(1,i);
			helices.push_back(temp_vec);

		}
		else if (this_val == ss3_HELIX){
			helices[num_helices-1].push_back(i);
			helix_length++;
		}


		last_val = this_val;
	}

	return helices;
}

int_vector_vector get_loops_vector(const PRODART::POSE::three_state_sec_struct_vector& secs){
	if (secs.size() == 0){
		return int_vector_vector();
	}

	int_vector_vector loops;

	//int_vector::const_iterator iter;
	loops.clear();

	three_state_sec_struct last_val = ss3_UNDEF;
	int num_loops = 0;
	int loop_length = 0;


	const int resNum = secs.size();

	for (int i = 0; i < resNum; i++){
		const three_state_sec_struct this_val = secs[i];

		if (last_val != ss3_OTHER && this_val == ss3_OTHER){
			num_loops++;
			loop_length = 1;
			int_vector temp_vec(1,i);
			loops.push_back(temp_vec);
			//PRINT_EXPR(i);
			//PRINT_EXPR(temp_vec.size());

		}
		else if (this_val == ss3_OTHER){
			loops[num_loops-1].push_back(i);
			loop_length++;
		}


		last_val = this_val;
	}

	return loops;
}


const double max_hb_dist = 5.5;
bool are_strands_paired(int_vector& strand1, int_vector& strand2, PRODART::POSE::pose_shared_ptr protein){
	int_vector::const_iterator iter_1, iter_2;
	int interact_count = 0;
	for (iter_1 = strand1.begin(); iter_1 != strand1.end(); iter_1++){
		const vector3d vec_1 = protein->get_bb_atom(POSE::CA, *iter_1)->get_coords(); //protein.getVecByAtomIndex(*iter_1);
		for (iter_2 = strand2.begin(); iter_2 != strand2.end(); iter_2++){

			const vector3d vec_2 = protein->get_bb_atom(POSE::CA, *iter_2)->get_coords(); //protein.getVecByAtomIndex(*iter_2);

			const double dist = (vec_1 - vec_2).mod();
			/*
			cout << (*iter_1 - PRODART::PROT::getRelAtomPosition(CA)) / PRODART::PROT::num_mainchain_atoms_per_residue << "\t"
				 << (*iter_2 - PRODART::PROT::getRelAtomPosition(CA)) / PRODART::PROT::num_mainchain_atoms_per_residue << "\t"
				 << dist << endl;
			*/
			if (dist < max_hb_dist){
				interact_count++;

			}


		}
	}

	if (interact_count >= 2){
		return true;
	}

	return false;
}

bool are_strands_parallel(int_vector& strand1, int_vector& strand2, PRODART::POSE::pose_shared_ptr protein){

	const int s1_start = strand1.front();
	const int s1_end = strand1.back();

	const int s2_start = strand2.front();
	const int s2_end = strand2.back();

	const vector3d s1_vec = protein->get_bb_atom(POSE::CA, s1_end)->get_coords() - protein->get_bb_atom(POSE::CA, s1_start)->get_coords();
	const vector3d s2_vec = protein->get_bb_atom(POSE::CA, s2_end)->get_coords() - protein->get_bb_atom(POSE::CA, s2_start)->get_coords();

	if ( s1_vec.dot(s2_vec) > 0 ){
		return true;
	}
	else {
		return false;
	}

}

const double max_strand_pairing_dist = 7.0;
residue_pair_map find_parallel_pairing(const int_vector& strand1, const int_vector& strand2,  PRODART::POSE::pose_shared_ptr protein){
	/*
	const int s1_start = strand1.front();
	const int s1_end = strand1.back();

	const int s2_start = strand2.front();
	const int s2_end = strand2.back();
	*/

	residue_pair_map rtn_map;

	const int s1_size = strand1.size();
	const int s2_size = strand2.size();

	const int max_neg = -(s1_size - 1);
	const int max_pos =  (s2_size - 1);

	int best_num_hb = 0, best_num_possible_hb = 0, best_relative_pos = 0;

	for (int i = max_neg; i <= max_pos; i++){
		int num_hb = 0;
		int num_possible_hb = 0;

		for (int index1 = 0; index1 < s1_size; index1++ ){
			int index2 = index1 + i;
			if (index1 >=0 && index2 >= 0 && index2 < s2_size){
				num_possible_hb++;
				const int at1 = strand1[index1];
				const int at2 = strand2[index2];
				const vector3d vec1 = protein->get_bb_atom(POSE::CA, at1)->get_coords();
				const vector3d vec2 = protein->get_bb_atom(POSE::CA, at2)->get_coords();

				const double dist = (vec1 - vec2).mod();

				if (dist < max_hb_dist){
					num_hb++;
				}


			}

		}

		if (num_hb > best_num_hb){
			best_relative_pos = i;
			best_num_hb = num_hb;
			best_num_possible_hb = num_possible_hb;
		}
		else if (num_hb == best_num_hb && num_possible_hb > best_num_possible_hb){
			best_relative_pos = i;
			best_num_hb = num_hb;
			best_num_possible_hb = num_possible_hb;
		}


	}

	int best_index1 = 0;
	if (best_relative_pos >= 0 ){
		best_index1 = 0;
	}
	else {
		best_index1 = - best_relative_pos;
	}
	const int best_index2 = best_index1 + best_relative_pos;

	//int_vec_iterator_pair return_pair;

	const int first_1 = strand1.front() + best_index1;
	const int first_2 = strand2.front() + best_index2;


	cout << "found parallel CA hbond pair:\t"
		 << first_1 << "\t"
		 << first_2 << "\t"
		 << best_relative_pos << "\t"
		 << best_num_hb << "\t"
		 << best_num_possible_hb << "\t"
		 << endl;


	int iter1 = first_1;
	int iter2 = first_2;
	while (iter1 <=  strand1.back() && iter2 <= strand2.back()){
		const int res1 = iter1;
		const int res2 = iter2;

		const vector3d vec1 = protein->get_bb_atom(POSE::CA, res1)->get_coords();
		const vector3d vec2 = protein->get_bb_atom(POSE::CA, res2)->get_coords();

		const double dist = (vec1 - vec2).mod();

		if (dist < max_strand_pairing_dist){

			rtn_map[protein->get_residue(res1)] = protein->get_residue(res2);

			/*
			cout << "CA hbond_restraint added:\t" << res1 << "\t"
				 << res2 << "\t"
				 << dist << "\t"
				 << endl;
				 */
		}
		else {
			/*
			cout << "CA hbond_restraint not added:\t" << res1  << "\t"
				 << res2 << "\t"
				 << dist << "\t"
				 << endl;
				 */
		}

		iter1++;
		iter2++;
	}

	return rtn_map;

}

residue_pair_map find_antiparallel_pairing(const int_vector& strand1, const int_vector& strand2,  PRODART::POSE::pose_shared_ptr protein){
	residue_pair_map rtn_map;


	/*
	const int s1_start = strand1.front();
	const int s1_end = strand1.back();

	const int s2_start = strand2.front();
	const int s2_end = strand2.back();
	*/

	const int s1_size = strand1.size();
	const int s2_size = strand2.size();

	const int max_neg = 0;//-(s1_size - 1);
	const int max_pos = (s1_size - 1) + (s2_size - 1);//(s2_size - 1);

	int best_num_hb = 0, best_num_possible_hb = 0, best_relative_pos = 0;

	for (int i = max_neg; i <= max_pos; i++){
		int num_hb = 0;
		int num_possible_hb = 0;

		for (int index1 = 0; index1 < s1_size; index1++ ){
			int index2 = i - index1;
			if (index1 >=0 && index2 >= 0 && index2 < s2_size){
				num_possible_hb++;
				const int res1 = strand1[index1];
				const int res2 = strand2[index2];
				const vector3d vec1 = protein->get_bb_atom(POSE::CA, res1)->get_coords();
				const vector3d vec2 = protein->get_bb_atom(POSE::CA, res2)->get_coords();

				const double dist = (vec1 - vec2).mod();

				if (dist < max_hb_dist){
					num_hb++;
				}


			}

		}

		if (num_hb > best_num_hb){
			best_relative_pos = i;
			best_num_hb = num_hb;
			best_num_possible_hb = num_possible_hb;
		}
		else if (num_hb == best_num_hb && num_possible_hb > best_num_possible_hb){
			best_relative_pos = i;
			best_num_hb = num_hb;
			best_num_possible_hb = num_possible_hb;
		}


	}


	int best_index2 = best_relative_pos;
	if (best_relative_pos < s2_size){
		best_index2 = best_relative_pos;
	}
	else {
		best_index2 = s2_size - 1;
	}
	const int best_index1 = best_relative_pos - best_index2;


	const int first_1 = strand1.front() + best_index1;
	const int first_2 = strand2.front() + best_index2;



	cout << "found anti-parallel CA hbond pair:\t"
		 << first_1 << "\t"
		 << first_2 << "\t"
		 << best_relative_pos << "\t"
		 << best_num_hb << "\t"
		 << best_num_possible_hb << "\t"
		 << endl;

	int iter1 = first_1;
	int iter2 = first_2;


	while (iter1 <= strand1.back() && iter2 >= strand2.front()){

		const int res1 = iter1;
		const int res2 = iter2;

		const vector3d vec1 = protein->get_bb_atom(POSE::CA, res1)->get_coords();
		const vector3d vec2 = protein->get_bb_atom(POSE::CA, res2)->get_coords();


		const double dist = (vec1 - vec2).mod();

		if (dist < max_strand_pairing_dist){
			rtn_map[protein->get_residue(res1)] = protein->get_residue(res2);

			/*
			cout << "CA hbond_restraint added:\t" << res1 << "\t"
				 << res2 << "\t"
				 << dist << "\t"
				 << endl;
				 */
		}
		else {
			/*
			cout << "CA hbond_restraint not added:\t" << res1 << "\t"
				 << res2 << "\t"
				 << dist << "\t"
				 << endl;
				 */
		}

		iter1++;
		iter2--;

	}

	return rtn_map;
}

bool add_strand_pair_rst(residue_pair_map& pairs, const bool is_parallel, const double weight){
	restraints_store* rst_store = restraints_store::Instance();
	residue_pair_map::iterator iter;

	for (iter = pairs.begin(); iter != pairs.end(); iter++){
		rst_store->add_strand_pair_rst(iter->first->get_trimmed_pdb_residue_index(), iter->first->get_chain()->getChainID(),
				iter->second->get_trimmed_pdb_residue_index(), iter->second->get_chain()->getChainID(),
				is_parallel, weight);
		cout << "STRAND_PAIR " << iter->first->get_trimmed_pdb_residue_index() << " " <<  iter->first->get_chain()->getChainID() << " "
				<< iter->second->get_trimmed_pdb_residue_index() << " " <<  iter->second->get_chain()->getChainID() << " ";
		if (is_parallel){
			cout << "PARALLEL ";
		}
		else {
			cout << "ANTIPARALLEL ";

		}
		const double dist = (iter->first->get_bb_atom(POSE::CA)->get_coords() - iter->second->get_bb_atom(POSE::CA)->get_coords()).mod();
		cout << weight << "\t# CA dist:\t" << dist << endl;
	}
	return true;
}

bool auto_add_strand_pair_restraints(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const double weight){
	if (secs.size() != (unsigned long)protein->get_residue_count()){
		cerr << "simple_ca2main_secs_refine_minimise: ERROR: sec struct (" << secs.size() << ") not same size as protein (" << protein->get_residue_count() << ")" << endl;
		return false;
	}

	int_vector_vector strands = get_strands_vector(secs);
	int_vector_vector helices = get_helices_vector(secs);
	const int num_strands = strands.size();
	const int num_helices = helices.size();
	if (num_strands == 0){
		cout << "\nsimple_ca2main_secs_refine_minimise: INFO: no beta strands defined\n" << endl;
		return true;
	}
	else {
		cout << "\nsimple_ca2main_secs_refine_minimise: INFO: " << num_strands
		<<" beta strands defined" << endl;
	}

	if (num_helices == 0){
		cout << "\nsimple_ca2main_secs_refine_minimise: INFO: no helices defined\n" << endl;
		//return true;
	}
	else {
		cout << "simple_ca2main_secs_refine_minimise: INFO: " << num_helices
		<<" alpha helices defined"
		<< endl
		<< endl;
	}

	int num_pairs = 0;
	int num_paral_pairs = 0;
	for (int strand_i = 0; strand_i < num_strands; strand_i++){
		for (int strand_j = strand_i + 1; strand_j < num_strands; strand_j++){

			if (are_strands_paired(strands[strand_i], strands[strand_j], protein)){
				num_pairs++;
				const bool is_parallel = are_strands_parallel(strands[strand_i], strands[strand_j], protein);
				if (is_parallel) {
					num_paral_pairs++;
					residue_pair_map pairing = find_parallel_pairing(strands[strand_i], strands[strand_j], protein);
					add_strand_pair_rst(pairing, is_parallel, weight);
				}
				else {
					residue_pair_map pairing = find_antiparallel_pairing(strands[strand_i], strands[strand_j], protein);
					add_strand_pair_rst(pairing, is_parallel, weight);
				}
			}

		}
	}
	cout << "\nsimple_ca2main_secs_refine_minimise: INFO: gaussian_distance_conserve: " << num_pairs
		 << " beta strand pairs found" << endl;
	cout << "simple_ca2main_secs_refine_minimise: INFO: of which " << num_paral_pairs
		 << " parallel strand pairs" << endl;
	cout << "simple_ca2main_secs_refine_minimise: INFO: of which " << (num_pairs - num_paral_pairs)
		 << " anti-parallel strand pairs" << endl;

	return true;
}

bool auto_add_secs_restraints(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const double weight){
	restraints_store* rst_store = restraints_store::Instance();
	const int resNum = protein->get_residue_count();
	for (int i = 0; i < resNum; i++){
		const three_state_sec_struct this_val = secs[i];
		if (this_val == ss3_HELIX || this_val == ss3_STRAND){
			residue_shared_ptr res = protein->get_residue(i);
			rst_store->add_sec_struct_rst(res->get_trimmed_pdb_residue_index(), res->get_chain()->getChainID(),
					to_4state(this_val), weight);
			cout << "SEC_STRUCT " << res->get_trimmed_pdb_residue_index() << " "
					<< res->get_chain()->getChainID() << " "
					<< to_char(to_4state(this_val)) << " " << weight
					<< endl;
		}
	}
	return true;
}

bool auto_add_secs_axis_restraints(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const double alpha_weight,
		const double beta_weight){
	restraints_store* rst_store = restraints_store::Instance();
	//const int resNum = protein->get_residue_count();
	int_vector_vector strands = get_strands_vector(secs);
	int_vector_vector helices = get_helices_vector(secs);
	const int num_strands = strands.size();
	const int num_helices = helices.size();

	for (int i = 0; i < num_strands; i++){
		const int start_num = strands[i].front();
		const int end_num = strands[i].back();
		//PRINT_EXPR(start_num);
		//PRINT_EXPR(end_num);

		if (end_num - start_num + 1 >= 4){
			vector3d start, end;
			get_strand_end_points_long(protein, start_num, end_num, start, end);
			vector3d noise(ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01));
			vector3d noise2(ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01));
			rst_store->add_sse_axis_rst(protein->get_residue(start_num)->get_trimmed_pdb_residue_index(), protein->get_residue(start_num)->get_chain()->getChainID(),
					protein->get_residue(end_num)->get_trimmed_pdb_residue_index(), protein->get_residue(end_num)->get_chain()->getChainID(),
					ss4_STRAND, beta_weight,
					start+ noise, end+noise2);
			cout << "SSE_AXIS STRAND " << protein->get_residue(start_num)->get_trimmed_pdb_residue_index() << " " << protein->get_residue(start_num)->get_chain()->getChainID() << " "
					<< protein->get_residue(end_num)->get_trimmed_pdb_residue_index() << " " << protein->get_residue(end_num)->get_chain()->getChainID() << " "
					<< beta_weight
					<< endl;
		}

	}

	for (int i = 0; i < num_helices; i++){
		const int start_num = helices[i].front();
		const int end_num = helices[i].back();
		//PRINT_EXPR(start_num);
		//PRINT_EXPR(end_num);

		if (end_num - start_num + 1 >= 4){
			vector3d start, end;
			get_helix_end_points_long(protein, start_num, end_num, start, end);
			vector3d noise(ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01));
			vector3d noise2(ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01), ENV::get_random_num_gen()->rand(0.01));
			rst_store->add_sse_axis_rst(protein->get_residue(start_num)->get_trimmed_pdb_residue_index(), protein->get_residue(start_num)->get_chain()->getChainID(),
					protein->get_residue(end_num)->get_trimmed_pdb_residue_index(), protein->get_residue(end_num)->get_chain()->getChainID(),
					ss4_HELIX, alpha_weight,
					start + noise, end + noise2);
			cout << "SSE_AXIS HELIX " << protein->get_residue(start_num)->get_trimmed_pdb_residue_index() << " " << protein->get_residue(start_num)->get_chain()->getChainID() << " "
					<< protein->get_residue(end_num)->get_trimmed_pdb_residue_index() << " " << protein->get_residue(end_num)->get_chain()->getChainID() << " "
					<< alpha_weight
					<< endl;
		}

	}

	return true;
}

bool auto_add_coord_restraints(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const double alpha_weight,
		const double beta_weight){
	restraints_store* rst_store = restraints_store::Instance();

	const int resNum = protein->get_residue_count();
	for (int i = 0; i < resNum; i++){
		const three_state_sec_struct this_val = secs[i];
		if (this_val == ss3_HELIX || this_val == ss3_STRAND){
			const double weight = this_val == ss3_HELIX ? alpha_weight : beta_weight;
			residue_shared_ptr res = protein->get_residue(i);
			const vector3d vec = res->get_bb_atom(POSE::CA)->get_coords();
			rst_store->add_coord_rst(res->get_trimmed_pdb_residue_index(), res->get_chain()->getChainID(), POSE::atom_type("CA"),
					vec, weight, 0.0);
			/*
			add_coord_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
					const UTILS::vector3d pos, const double weight, const double equil_dist = 0);
			 */

			cout << "COORD_HARMONIC " << res->get_trimmed_pdb_residue_index() << " "
					<< res->get_chain()->getChainID() << " "
					<< "CA "
					<< vec.x << " "
					<< vec.y << " "
					<< vec.z << " "
					<< weight << " "
					<< 0.0
					<< endl;
		}
	}
	return true;


}

bool simple_ca2main_secs_refine_minimise(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const unsigned long ca_ref_steps,
		const unsigned long min_steps,
		const std::string bb_min_potentials_weight_set){

	if (secs.size() != (unsigned long)protein->get_residue_count()){
		cerr << "simple_ca2main_secs_refine_minimise: ERROR: sec struct (" << secs.size() << ") not same size as protein (" << protein->get_residue_count() << ")" << endl;
		return false;
	}


	bool result = auto_add_strand_pair_restraints(protein, secs, 100.0);
	if (result) result = auto_add_secs_restraints(protein, secs, 10.0);
	if (result) result = auto_add_coord_restraints(protein, secs, 0.5*0.1, 0.5*0.025);
	if (result) result = auto_add_secs_axis_restraints(protein, secs, 1, 0.1);
	if (result) result = LOOP::ca_minimise_bonds_bumps(protein);
	if (result) result = ca_mc_refine(protein, ca_ref_steps, 1.0, "ca_default");
	if (result) result = new_ca2main_minimise(protein, min_steps / 10);
	cout << endl;
	if (result) result = bb_minimise(protein, "bb_min_no_strand_hb_only_sp_rst", min_steps / 5);
	cout << endl;
	for (int i = 0; i < 5; i++){
		if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps / 10);
		cout << endl;
		if (result) result = bb_minimise(protein, "bb_min_no_strand_hb_only_sp_rst", min_steps / 5);
		cout << endl;
		if (result) result = bb_minimise(protein, "bb_min_default", min_steps / 10);
		cout << endl;
	}
	if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps);
	cout << endl;
	/*
	if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps / 50);
	cout << endl;
	if (result) result = bb_minimise(protein, "bb_min_default", min_steps / 50);
	cout << endl;
	if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps / 10);
	cout << endl;
	if (result) result = bb_minimise(protein, "bb_min_default", min_steps);
	 */
	//cerr << "simple_ca2main_secs_refine_minimise: ERROR: not implemented yet" << endl;
	return result;

}

bool ca_mc_refine(PRODART::POSE::pose_shared_ptr protein, const unsigned long steps, const double beta, const std::string ca_pot){

	potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot = pot_factory->make_preset_potentials_container(ca_pot);

	PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(protein);

	PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = MOVERS::move_set_factory::Instance()->make_preset_move_set(meta_data,
			"ca_standard",
			PRODART::ENV::get_random_num_gen());


	PRODART::POSE::SIM::monte_carlo_shared_ptr mc_sim = PRODART::POSE::SIM::new_monte_carlo(pot,
			overall_move_set,
			PRODART::POSE::SIM::CA::new_simple_ca_mc_protocol(steps,
					beta,
					PRODART::ENV::get_random_num_gen()));


	cout << "PROTOCOLS::CA2MAIN::ca_mc_refine: running ca mc for " << steps << " steps..." << endl;

	const bool result = mc_sim->make_move(meta_data);
	mc_sim->get_state().print_summary(cout);
	cout << "done" << endl;
	return result;
}

bool activate_sse_only(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const int margin){

	const int res_count = protein->get_residue_count();

	int act_count = 0;
	int inact_count =0;
	for (int i = 0; i < res_count; i++){
		residue_shared_ptr res = protein->get_residue(i);
		const three_state_sec_struct s_i = secs[i];
		const three_state_sec_struct s_i_m_m = i - margin > 0 ? secs[i - margin] : ss3_OTHER;
		const three_state_sec_struct s_i_p_m = i + margin < res_count ? secs[i + margin] : ss3_OTHER;

		const bool activate = (s_i != ss3_OTHER || s_i_m_m != ss3_OTHER || s_i_p_m != ss3_OTHER);

		atom_shared_ptr_vector all_atoms = res->get_all_atoms();
		for (unsigned int j = 0; j < all_atoms.size(); j++){
			all_atoms[j]->setActive(activate && all_atoms[j]->isActive());
			if (activate && all_atoms[j]->isActive()){
				act_count++;
			}
			else {
				inact_count++;
			}
		}

	}
	PRINT_EXPR(act_count);
	PRINT_EXPR(inact_count);

	return true;
}

PRODART::POSE::residue_shared_ptr_vector add_dummy_loop_residues(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const int num_res,
		PRODART::POSE::three_state_sec_struct_vector& mod_secs,
		const double end_rst_wt){

	restraints_store* rst_store = restraints_store::Instance();

	mod_secs = secs;
	int_vector_vector loops = get_loops_vector(secs);
	PRODART::POSE::residue_shared_ptr_vector dummyres;
	PRODART::POSE::residue_type_vector rtVec(num_res, residue_type("GLY"));
	const PRODART::POSE::three_state_sec_struct_vector insert_secs(num_res, ss3_OTHER);
	for (int i = loops.size()-1; i >=0; i--){
		const int start_num = loops[i].front();
		const int end_num = loops[i].back();
		const int len = loops[i].size();

		if (start_num != 0 && end_num != protein->get_residue_count()-1){
			const int insert_res = len / 2 + start_num;
			residue_shared_ptr res1 = protein->get_residue(insert_res);
			residue_shared_ptr res2 = protein->get_residue(insert_res-1);
			rst_store->add_bond_harm_rst(res1->get_pdb_residue_index(), res1->get_chain()->getChainID(), atom_type("CA"),
					res2->get_pdb_residue_index(), res2->get_chain()->getChainID(), atom_type("CA"),
					3.8, end_rst_wt);
			PRINT_EXPR(insert_res);
			PRODART::POSE::residue_shared_ptr_vector temp =LOOP::insert_ca_loop_linear_interpolate(protein, insert_res,rtVec);
			dummyres.insert(dummyres.end(), temp.begin(), temp.end());
			mod_secs.insert(mod_secs.begin()+insert_res, insert_secs.begin(), insert_secs.end());
			for (unsigned int j = 0; j < temp.size(); j++){
				temp[j]->get_bb_atom(POSE::CA)->setActive(false);
			}
		}

	}

	protein->index();
	return dummyres;
}

bool delete_dummy_loop_residues(PRODART::POSE::pose_shared_ptr protein,
		PRODART::POSE::residue_shared_ptr_vector res_to_delete){
	for (unsigned int i = 0; i < res_to_delete.size(); i++){
		protein->delete_residue(res_to_delete[i]->get_internal_residue_index());
	}
	return true;
}

// simple ca2main with sec struct refinement and minimisation of loops and sec struct independently
bool simple_ca2main_secs_indep_refine_minimise(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::three_state_sec_struct_vector& secs,
		const unsigned long min_steps,
		const std::string bb_min_potentials_weight_set){

	cerr << "\nsimple_ca2main_secs_indep_refine_minimise: ERROR: THIS IS BUGGY DO NOT USE YET!!!\n" << endl;


	if (secs.size() != (unsigned long)protein->get_residue_count()){
		cerr << "simple_ca2main_secs_refine_minimise: ERROR: sec struct (" << secs.size() << ") not same size as protein (" << protein->get_residue_count() << ")" << endl;
		return false;
	}


	const int min_cycles = 7;

	PRODART::POSE::three_state_sec_struct_vector mod_secs;

	const LOOP::atom_shared_ptr_bool_map orig_act_map = LOOP::get_activity_map(protein);
	//activate_sse_only(protein, secs,4);

	PRINT_EXPR(protein->get_residue_count());
	PRODART::POSE::residue_shared_ptr_vector  dummy_loop_res = add_dummy_loop_residues(protein, secs, 5, mod_secs, 0.5);
	PRINT_EXPR(protein->get_residue_count());
	PRINT_EXPR(secs);
	PRINT_EXPR(mod_secs);

	bool result = auto_add_strand_pair_restraints(protein, mod_secs, 100.0);
	if (result) result = auto_add_secs_restraints(protein, mod_secs, 10.0);
	if (result) result = auto_add_coord_restraints(protein, mod_secs, 0.5*0.1, 0.5*0.010);
	if (result) result = auto_add_secs_axis_restraints(protein, mod_secs, 0.5, 0.05);
	if (result) result = LOOP::ca_minimise_bonds_bumps(protein);
	if (result) result = ca_mc_refine(protein, 100000, 1.0, "ca_default");
	if (result) result = new_ca2main_minimise(protein, min_steps / 10);
	cout << endl;
	for (int i = 0; i < min_cycles; i++){
		if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps / 10);
		cout << endl;
		if (result) result = bb_minimise(protein, "bb_min_default", min_steps / 10);
		cout << endl;
	}
	if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps);
	cout << endl;

	delete_dummy_loop_residues(protein, dummy_loop_res);
	LOOP::restore_activity_map(orig_act_map);

	bool_vector loop_mask = make_residue_loop_mask(secs);

	LOOP::run_ca_loop_anneal(protein, loop_mask, "ca_default", 50, 10000, 100000, 0.1, 1.2);
	new_ca2main_minimise(protein, loop_mask, 500, "bb_min_default");

	for (int i = 0; i < min_cycles; i++){
		if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps / 10);
		cout << endl;
		if (result) result = bb_minimise(protein, "bb_min_default", min_steps / 10);
		cout << endl;
	}
	if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps);
	cout << endl;

	/*
	if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps / 50);
	cout << endl;
	if (result) result = bb_minimise(protein, "bb_min_default", min_steps / 50);
	cout << endl;
	if (result) result = bb_minimise(protein, "bb_min_default_no_rst", min_steps / 10);
	cout << endl;
	if (result) result = bb_minimise(protein, "bb_min_default", min_steps);
	 */
	//cerr << "simple_ca2main_secs_refine_minimise: ERROR: not implemented yet" << endl;
	cerr << "\nsimple_ca2main_secs_indep_refine_minimise: ERROR: THIS IS BUGGY DO NOT USE YET!!!\n" << endl;
	return result;
}


// ca2main with no CA refinement, minimisation or optimisation with new alphabet method
bool alphabet_ca2main_fixed_ca_minimise(PRODART::POSE::pose_shared_ptr protein,
		const PRODART::POSE::BB_BUILDER::alphabet_bb_builder_shared_ptr bbb,
		const unsigned long min_steps,
		const std::string bb_min_potentials_weight_set){

	const double res_b_fact = 1.0;

	for (int i = 0 ; i < protein->get_all_atom_count(); i++ ){
		atom_shared_ptr atm = protein->get_atom(i);
		atm->set_b_factor(res_b_fact);
	}

	PRODART::POSE::pose_shared_ptr orig_pose = protein->clone();

	/* // outputs each fitter fragments
	for (int i = 0; i < protein->get_residue_count()-5; i++){
		pose_shared_ptr letter = bbb->get_best_fit_letter(i, protein);
		//bbb->build_fragment(i, protein);

		letter->outputPdb(cout);
	}
	*/

	CA2MAIN::add_two_end_ca_atoms_residues(protein);
	CA2MAIN::add_two_end_ca_atoms_residues(protein);

	atom_shared_ptr_vector added_atoms = bbb->build_backbone(protein);

	CA2MAIN::delete_two_end_residues(protein);
	CA2MAIN::delete_two_end_residues(protein);

	std::vector<boost::weak_ptr<POSE::atom> > weak_added_atoms;

	for (atom_shared_ptr_vector::iterator iter = added_atoms.begin(); iter != added_atoms.end(); iter++){
		weak_added_atoms.push_back(*iter);
	}
	added_atoms.clear();

	for (std::vector<boost::weak_ptr<POSE::atom> >::iterator iter = weak_added_atoms.begin(); iter != weak_added_atoms.end(); iter++){
		if (!(*iter).expired()){
			added_atoms.push_back((*iter).lock());
		}
	}

	//minimise here
	if (min_steps !=0){

		bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);

		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container(bb_min_potentials_weight_set);
		potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();
		bb_min_pot->get_energy(bb_meta,enrg_map);
		enrg_map.print_headers(cout);
		enrg_map.print_weighted_components(cout);

		bool_vector atom_selection(protein->get_all_atom_count(), false);
		for (atom_shared_ptr_vector::const_iterator iter = added_atoms.begin(); iter != added_atoms.end(); iter++){
			if ((*iter)->isActiveAndSet()){
				atom_selection[(*iter)->get_seq_num()] = true;
			}
		}

		minimiser_interface_shared_ptr min_sim =  new_minimiser(bb_min_pot, atom_selection);
		min_sim->set_steps(min_steps);

		//const bool result =
		min_sim->make_move(bb_meta);

		bb_min_pot->get_energy(bb_meta,enrg_map);
		enrg_map.print_headers(cout);
		enrg_map.print_weighted_components(cout);

	}



	//



	const prot_backbone_map* bb_map = POSE::prot_backbone_map::Instance();

	const std::vector<BBAtomType> bb_vec =  bb_map->get_smaller_bb_atom_list();


	for (int i = 0; i < protein->get_residue_count(); i++){
		residue_shared_ptr res_rb = protein->get_residue(i);

		if (protein->get_bb_atom(POSE::N, i)->get_b_factor() != res_b_fact
				|| protein->get_bb_atom(POSE::C, i)->get_b_factor() != res_b_fact){
			for (std::vector<BBAtomType>::const_iterator iter = bb_vec.begin() ; iter != bb_vec.end(); iter++){
				atom_shared_ptr nat = protein->get_bb_atom(*iter, i);
				atom_shared_ptr origat = orig_pose->get_bb_atom(*iter, i);
				if (nat && origat){
					if (nat->get_b_factor() != res_b_fact
							&& nat->isActiveAndSet() && origat->isActiveAndSet()){
						stringstream appline;
						appline << "ATOM_DIST\t" << res_rb->get_trimmed_pdb_residue_index() << "\t"
								<< res_rb->get_chain()->getChainID() << "\t"
								<< nat->get_type().get_trimmed_label() << "\t"
								<< (nat->get_coords() - origat->get_coords()).mod() << "\t"
								<< res_rb->get_type().get_label3() << "\t"
								<< endl;
						protein->add_appendline(appline.str());
					}
				}
			}
			const double new_phi =  UTILS::radians_to_degrees(protein->get_phi(i));
			const double new_psi =   UTILS::radians_to_degrees(protein->get_psi(i));
			const double new_omega =  UTILS::radians_to_degrees(protein->get_omega_to_prev(i));
			residue_shared_ptr res_orig = orig_pose->get_residue(i);
			const double orig_phi =  UTILS::radians_to_degrees(orig_pose->get_phi(i));
			const double orig_psi =   UTILS::radians_to_degrees(orig_pose->get_psi(i));
			const double orig_omega =  UTILS::radians_to_degrees(orig_pose->get_omega_to_prev(i));
			stringstream appline;
			appline << "DIHEDRALS\t" << res_rb->get_trimmed_pdb_residue_index() << "\t"
					<< res_rb->get_chain()->getChainID() << "\t"
					<< new_phi << "\t"
					<< new_psi << "\t"
					<< new_omega << "\t"
					<< orig_phi << "\t"
					<< orig_psi << "\t"
					<< orig_omega << "\t"
					<< res_rb->get_type().get_label3() << "\t"
					<< endl;
			protein->add_appendline(appline.str());
		}

	}
	return true;
}


double ca2main_get_rmsd_special(PRODART::POSE::pose_shared_ptr ref, PRODART::POSE::pose_shared_ptr protein){
	const int resCount1 = ref->get_residue_count();
	const int resCount2 = protein->get_residue_count();
	if (resCount1 != resCount2){
		std::cerr << "ERROR: ca2main_get_rmsd_special can not calculate RMSD between different size proteins " <<  protein->get_label() << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	const int chainCount1 = ref->get_chain_count();
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "CHAIN_COUNT\t" << chainCount1 << endl;
	const int chainCount2 = protein->get_chain_count();
	if (chainCount1 != chainCount2){
		std::cerr << "ERROR: ca2main_get_rmsd_special can not calculate RMSD between different chain counts" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	const double ca_rmsd = PRODART::POSE_UTILS::get_ca_rmsd_superpose(ref, protein);
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "CA_RMSD\t" << ca_rmsd << "\t" << resCount1 << endl;

	vector<atom_shared_ptr> n_vec_ref(0), n_vec_prot(0);
	vector<atom_shared_ptr> c_vec_ref(0), c_vec_prot(0);
	vector<atom_shared_ptr> o_vec_ref(0), o_vec_prot(0);
	vector<atom_shared_ptr> cb_vec_ref(0), cb_vec_prot(0);

	vector<atom_shared_ptr> all_vec_ref(0), all_vec_prot(0);

	MISC::clean_pose(protein);

	for (int ch = 0; ch < chainCount1; ch++ ){
		chain_shared_ptr ref_ch = ref->get_chain(ch);
		chain_shared_ptr prot_ch = protein->get_chain(ch);
		std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping;
		PRODART::POSE::atom_shared_ptr_vector atoms_to_move(0);
		if (ref_ch->length() != prot_ch->length()){
			std::cerr << "ERROR: ca2main_get_rmsd_special chain lengths are different size" << endl;
			return std::numeric_limits<double>::quiet_NaN();
		}

		for (int i = 0; i < ref_ch->length(); i++){

			atom_shared_ptr n_r_at = ref_ch->get_bb_atom(N,i);
			atom_shared_ptr c_r_at = ref_ch->get_bb_atom(C,i);
			atom_shared_ptr o_r_at = ref_ch->get_bb_atom(O,i);
			atom_shared_ptr cb_r_at = ref_ch->get_bb_atom(CB,i);
			atom_shared_ptr ca_r_at = ref_ch->get_bb_atom(POSE::CA,i);

			atom_shared_ptr n_p_at = prot_ch->get_bb_atom(N,i);
			atom_shared_ptr c_p_at = prot_ch->get_bb_atom(C,i);
			atom_shared_ptr o_p_at = prot_ch->get_bb_atom(O,i);
			atom_shared_ptr cb_p_at = prot_ch->get_bb_atom(CB,i);
			atom_shared_ptr ca_p_at = prot_ch->get_bb_atom(POSE::CA,i);


			atom_mapping[ca_r_at] = ca_p_at;

			atom_shared_ptr_vector all_ats = prot_ch->get_residue(i)->get_all_atoms();
			atoms_to_move.insert( atoms_to_move.end(), all_ats.begin(), all_ats.end() );

			double ref_phi = 0, ref_psi = 0, ref_omega = 0;
			double phi = 0, psi = 0, omega = 0;
			bool dihs_set = false;

			if (i == 0){
				// nothing
			}
			else if (i == 1){
				// nothing
				c_vec_ref.push_back(c_r_at);
				c_vec_prot.push_back(c_p_at);
				o_vec_ref.push_back(o_r_at);
				o_vec_prot.push_back(o_p_at);

				ref_phi = ref_ch->get_residue(i)->get_phi();
				ref_psi = ref_ch->get_residue(i)->get_psi();
				ref_omega = ref_ch->get_residue(i)->get_omega_to_prev();
				phi = prot_ch->get_residue(i)->get_phi();
				psi = prot_ch->get_residue(i)->get_psi();
				omega = prot_ch->get_residue(i)->get_omega_to_prev();
				dihs_set = true;
			}
			else if (i == ref_ch->length() -1 ){
				// nothing
			}
			else if (i == ref_ch->length() -2 ){
				// nothing
				c_vec_ref.push_back(c_r_at);
				c_vec_prot.push_back(c_p_at);
				o_vec_ref.push_back(o_r_at);
				o_vec_prot.push_back(o_p_at);

				ref_phi = ref_ch->get_residue(i)->get_phi();
				ref_psi = ref_ch->get_residue(i)->get_psi();
				ref_omega = ref_ch->get_residue(i)->get_omega_to_prev();
				phi = prot_ch->get_residue(i)->get_phi();
				psi = prot_ch->get_residue(i)->get_psi();
				omega = prot_ch->get_residue(i)->get_omega_to_prev();
				dihs_set = true;
			}
			else if (i == 2){
				c_vec_ref.push_back(c_r_at);
				c_vec_prot.push_back(c_p_at);
				o_vec_ref.push_back(o_r_at);
				o_vec_prot.push_back(o_p_at);

				ref_phi = ref_ch->get_residue(i)->get_phi();
				ref_psi = ref_ch->get_residue(i)->get_psi();
				ref_omega = ref_ch->get_residue(i)->get_omega_to_prev();
				phi = prot_ch->get_residue(i)->get_phi();
				psi = prot_ch->get_residue(i)->get_psi();
				omega = prot_ch->get_residue(i)->get_omega_to_prev();
				dihs_set = true;

			}
			else if (i == ref_ch->length() - 3){
				n_vec_ref.push_back(n_r_at);
				n_vec_prot.push_back(n_p_at);

				ref_phi = ref_ch->get_residue(i)->get_phi();
				ref_psi = ref_ch->get_residue(i)->get_psi();
				ref_omega = ref_ch->get_residue(i)->get_omega_to_prev();
				phi = prot_ch->get_residue(i)->get_phi();
				psi = prot_ch->get_residue(i)->get_psi();
				omega = prot_ch->get_residue(i)->get_omega_to_prev();
				dihs_set = true;
			}
			else {
				n_vec_ref.push_back(n_r_at);
				n_vec_prot.push_back(n_p_at);
				c_vec_ref.push_back(c_r_at);
				c_vec_prot.push_back(c_p_at);
				o_vec_ref.push_back(o_r_at);
				o_vec_prot.push_back(o_p_at);
				if (cb_r_at && cb_p_at){
					if (cb_r_at->isActiveAndSet() && cb_p_at->isActiveAndSet()){
						cb_vec_ref.push_back(cb_r_at);
						cb_vec_prot.push_back(cb_p_at);
					}
				}
				ref_phi = ref_ch->get_residue(i)->get_phi();
				ref_psi = ref_ch->get_residue(i)->get_psi();
				ref_omega = ref_ch->get_residue(i)->get_omega_to_prev();
				phi = prot_ch->get_residue(i)->get_phi();
				psi = prot_ch->get_residue(i)->get_psi();
				omega = prot_ch->get_residue(i)->get_omega_to_prev();
				dihs_set = true;
			}

			if (dihs_set){

				bool is_trans = (POSE_UTILS::get_phi_psi_omega_sector(ref_phi, ref_psi, ref_omega) !=  POSE::ss4_CIS);
// (POSE_UTILS::get_phi_psi_omega_sector(phi, psi, omega) !=  POSE::ss4_CIS)

				string transStr = is_trans ? string("TRANS") : string("CIS");

				const double dist = sqrt(pow(UTILS::get_dihedral_distance(phi, ref_phi),2) + pow(UTILS::get_dihedral_distance(psi, ref_psi),2));

				string distStr = UTILS::radians_to_degrees(dist) > 45 ? string("BIG_DIFF") : string("SMALL_DIFF");

				const double omega_dist = UTILS::get_dihedral_distance(omega, ref_omega);
				string omega_distStr = UTILS::radians_to_degrees(omega_dist) > 90 ? string("OMEGA_FLIP") : string("OMEGA_OK");


				cout  << fixed << setprecision(1)
						<< protein->get_label() << "\t"
						<< ref_ch->get_residue(i)->get_type().get_label3() << "\t"
						<< ref_ch->get_residue(i)->get_trimmed_pdb_residue_index() << "\t"
						<< ref_ch->getChainID() << "\t"
						<< "DIHEDRALS\t"
						<< transStr << "\t"
						<< UTILS::radians_to_degrees(ref_phi) << "\t"
						<< UTILS::radians_to_degrees(ref_psi) << "\t"
						<< UTILS::radians_to_degrees(ref_omega) << "\t"
						<< "rebuilt:\t"
						<< UTILS::radians_to_degrees(phi) << "\t"
						<< UTILS::radians_to_degrees(psi) << "\t"
						<< UTILS::radians_to_degrees(omega) << "\t"
						<< "2D_RAMA_DIST:\t" << UTILS::radians_to_degrees(dist) << "\t"
						<< distStr << "\t"
						<< "OMEGA_DIST:\t" << UTILS::radians_to_degrees(omega_dist) << "\t"
						<< omega_distStr << "\t"
						<< endl;
			}

		}
		PRODART::POSE_UTILS::get_rmsd_superpose(atoms_to_move, atom_mapping);
		//PRODART::POSE_UTILS::get_ca_rmsd_superpose(ref, protein, residue_mapping);

	}

	const double ca_again_rmsd = PRODART::POSE_UTILS::get_ca_rmsd_superpose(ref, protein);
	cout << protein->get_label() << "\t" << "CA_again_RMSD\t" << ca_again_rmsd << "\t" << resCount1 << endl;

	all_vec_ref.insert( all_vec_ref.end(), n_vec_ref.begin(), n_vec_ref.end() );
	all_vec_ref.insert( all_vec_ref.end(), c_vec_ref.begin(), c_vec_ref.end() );
	all_vec_ref.insert( all_vec_ref.end(), o_vec_ref.begin(), o_vec_ref.end() );

	all_vec_prot.insert( all_vec_prot.end(), n_vec_prot.begin(), n_vec_prot.end() );
	all_vec_prot.insert( all_vec_prot.end(), c_vec_prot.begin(), c_vec_prot.end() );
	all_vec_prot.insert( all_vec_prot.end(), o_vec_prot.begin(), o_vec_prot.end() );

	const double n_rmsd = get_fixed_position_rmsd(n_vec_ref, n_vec_prot);
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "N_RMSD\t" << n_rmsd << "\t" << n_vec_prot.size() << endl;
	const double c_rmsd = get_fixed_position_rmsd(c_vec_ref, c_vec_prot);
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "C_RMSD\t" << c_rmsd << "\t" << c_vec_prot.size() << endl;
	const double o_rmsd = get_fixed_position_rmsd(o_vec_ref, o_vec_prot);
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "O_RMSD\t" << o_rmsd << "\t" << o_vec_prot.size() << endl;
	const double cb_rmsd = get_fixed_position_rmsd(cb_vec_ref, cb_vec_prot);
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "CB_RMSD\t" << cb_rmsd << "\t" << cb_vec_prot.size() << endl;

	const double bb_rmsd = get_fixed_position_rmsd(all_vec_ref, all_vec_prot);
	cout << setprecision(3) << protein->get_label() << "\t" << "OVERALL_NCO_RMSD\t" << bb_rmsd << "\t" << all_vec_prot.size() << endl;


	cout << endl;

	return bb_rmsd;
}

//!get RMSD for ca2main method - special for Ben's ca2main paper. sidechain only (not including CB)
double ca2main_get_rmsd_sidechain(PRODART::POSE::pose_shared_ptr ref, PRODART::POSE::pose_shared_ptr protein){
	const int resCount1 = ref->get_residue_count();
	const int resCount2 = protein->get_residue_count();
	if (resCount1 != resCount2){
		std::cerr << "ERROR: ca2main_get_rmsd_special can not calculate RMSD between different size proteins " <<  protein->get_label() << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	const int chainCount1 = ref->get_chain_count();
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "CHAIN_COUNT\t" << chainCount1 << endl;
	const int chainCount2 = protein->get_chain_count();
	if (chainCount1 != chainCount2){
		std::cerr << "ERROR: ca2main_get_rmsd_special can not calculate RMSD between different chain counts" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	const double ca_rmsd = PRODART::POSE_UTILS::get_ca_rmsd_superpose(ref, protein);
	cout  << fixed << setprecision(3) << protein->get_label() << "\t" << "CA_RMSD\t" << ca_rmsd << "\t" << resCount1 << endl;

	MISC::clean_pose(protein);


	vector<atom_shared_ptr> sc_vec_ref(0), sc_vec_prot(0);
	for (int i =0; i < resCount1; i++){
		residue_shared_ptr r_res = ref->get_residue(i);
		residue_shared_ptr p_res = protein->get_residue(i);
		sidechain_shared_ptr r_sc = r_res->get_sidechain();
		atom_shared_ptr_vector all_sc_ats = r_sc->get_all_sc_atoms();
		vector<atom_shared_ptr> this_sc_vec_ref(0), this_sc_vec_prot(0);
		for (unsigned int s = 0; s < all_sc_ats.size(); s++){
			atom_shared_ptr r_at = all_sc_ats[s];
			atom_shared_ptr p_at = p_res->get_atom(r_at->get_type());
			if (r_at && p_at){
				if (r_at->isActiveAndSet() && p_at->isActiveAndSet()){
					sc_vec_ref.push_back(r_at);
					sc_vec_prot.push_back(p_at);

					this_sc_vec_ref.push_back(r_at);
					this_sc_vec_prot.push_back(p_at);
				}
			}
			else {
				cerr << "WARNING: atom not found in decoy: " << protein->get_label() << " " << r_at->get_type().get_trimmed_label() << endl;
			}
		}

		if (this_sc_vec_ref.size() != 0 && this_sc_vec_prot.size() != 0 ){
			atom_shared_ptr n_r_at = ref->get_bb_atom(N,i);
			atom_shared_ptr c_r_at = ref->get_bb_atom(C,i);
			atom_shared_ptr o_r_at = ref->get_bb_atom(O,i);
			atom_shared_ptr cb_r_at = ref->get_bb_atom(CB,i);
			atom_shared_ptr ca_r_at = ref->get_bb_atom(POSE::CA,i);

			atom_shared_ptr n_p_at = protein->get_bb_atom(N,i);
			atom_shared_ptr c_p_at = protein->get_bb_atom(C,i);
			atom_shared_ptr o_p_at = protein->get_bb_atom(O,i);
			atom_shared_ptr cb_p_at = protein->get_bb_atom(CB,i);
			atom_shared_ptr ca_p_at = protein->get_bb_atom(POSE::CA,i);

			vector<atom_shared_ptr> this_nco_vec_ref(0), this_nco_vec_prot(0);

			this_nco_vec_ref.push_back(n_r_at);
			this_nco_vec_ref.push_back(c_r_at);
			this_nco_vec_ref.push_back(o_r_at);

			this_nco_vec_prot.push_back(n_p_at);
			this_nco_vec_prot.push_back(c_p_at);
			this_nco_vec_prot.push_back(o_p_at);

			const double this_nco_rmsd = get_fixed_position_rmsd(this_nco_vec_ref, this_nco_vec_prot);

			const double this_sc_rmsd = get_fixed_position_rmsd(this_sc_vec_ref, this_sc_vec_prot);

			cout << setprecision(3) << protein->get_label() << "\t" << protein->get_residue(i)->get_type().get_label3() << "\t" << "per_res_NCO_SC_RMSD\t" << this_nco_rmsd << "\t" << this_sc_rmsd  << endl;

		}


	}

	MISC::clean_pose(protein);

	const double sc_rmsd = get_fixed_position_rmsd(sc_vec_ref, sc_vec_prot);
	cout << setprecision(3) << protein->get_label() << "\t" << "OVERALL_SC_RMSD\t" << sc_rmsd << "\t" << sc_vec_prot.size() << endl;

	return sc_rmsd;
}

}
}
}


