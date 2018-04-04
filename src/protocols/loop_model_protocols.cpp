/*
 * loop_model_protocols.cpp
 *
 *  Created on: 28 Sep 2010
 *      Author: jmacdona
 */
#include "loop_model_protocols.h"
#include "protocols.h"


using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE::SIM;
using namespace PRODART::POSE::POTENTIALS;

namespace PRODART {
namespace PROTOCOLS{
namespace LOOP{

PRODART::POSE::residue_shared_ptr_vector insert_blank_residues(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec){

	const int insert_size = rtVec.size();

	if (insert_size == 0){
		return residue_shared_ptr_vector(0);
	}

	//const int n_anchor_num = insert_resnum - 1;
	//const int c_anchor_num = insert_resnum;

	if (insert_resnum < 0 || insert_resnum > protein->get_residue_count()){
		return residue_shared_ptr_vector(0);
	}

	//residue_shared_ptr n_anchor_res = protein->get_residue(n_anchor_num);
	//residue_shared_ptr c_anchor_res = protein->get_residue(c_anchor_num);

	residue_shared_ptr_vector new_residues(insert_size);//= chainPtr->insertNewResidueRange(insert_resnum, rtVec);
	//protein.reindex();

	for (int i = 0; i < insert_size; i++){

		new_residues[i] = protein->insert_residue(rtVec[i],insert_resnum+i);
		new_residues[i]->set_pdb_residue_index("XXX");

	}

	protein->index();

	return new_residues;
}

PRODART::POSE::residue_shared_ptr_vector insert_ca_loop_linear_interpolate(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const double random_component,
		const double CoM_component ){

	const double epsilon = std::numeric_limits<double>::min();

	protein->index();

	MTRand::MTRand_shared_ptr randomNumGen = PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id());

	const int insert_size = rtVec.size();

	if (insert_size == 0){
		return residue_shared_ptr_vector(0);
	}

	const int n_anchor_num = insert_resnum - 1;
	const int c_anchor_num = insert_resnum;



	if (insert_resnum < 0 || insert_resnum > protein->get_residue_count()){
		return residue_shared_ptr_vector(0);
	}

	residue_shared_ptr n_anchor_res = protein->get_residue(n_anchor_num);
	residue_shared_ptr c_anchor_res = protein->get_residue(c_anchor_num);

	atom_shared_ptr n_anchor_atom = protein->get_bb_atom(POSE::CA, n_anchor_num);
	atom_shared_ptr c_anchor_atom = protein->get_bb_atom(POSE::CA, c_anchor_num);

	const vector3d n_anchor_vec = n_anchor_atom->get_coords();
	const vector3d c_anchor_vec = c_anchor_atom->get_coords();

	//cout << n_anchor_vec << endl;
	//cout << c_anchor_vec << endl;

	const double line_len = (c_anchor_vec - n_anchor_vec).mod();
	const double bond_len = line_len / static_cast<double>(insert_size+1);
	const vector3d unit_vec = (c_anchor_vec - n_anchor_vec) / line_len;

	const vector3d CoM = PRODART::POSE_UTILS::get_ca_CoM(protein);

	//cout << line_len << endl;
	//cout << bond_len << endl;
	//cout << unit_vec << endl;

	residue_shared_ptr_vector new_residues(insert_size);//= chainPtr->insertNewResidueRange(insert_resnum, rtVec);
	//protein.reindex();

	for (int i = 0; i < insert_size; i++){

		new_residues[i] = protein->insert_residue(rtVec[i],insert_resnum+i);

		vector3d new_ca_vec = n_anchor_vec + ((static_cast<double>(i+1) * bond_len) * unit_vec);
		//Atom* new_atom =
		//cout << i << "\t" << new_ca_vec << endl;
		const double vec_len = (new_ca_vec - CoM).mod();
		const vector3d from_com_vec = vec_len > epsilon ?(new_ca_vec - CoM) / (new_ca_vec - CoM).mod() : vector3d(0,0,0);
		/*
		PRINT_EXPR((new_ca_vec - CoM));
		PRINT_EXPR((new_ca_vec - CoM));
		PRINT_EXPR((new_ca_vec - CoM).mod());
		PRINT_EXPR((new_ca_vec - CoM) / (new_ca_vec - CoM).mod());
		PRINT_EXPR(from_com_vec);
		*/

		if (random_component != 0){
			const vector3d transVec(randomNumGen->randNorm(0, random_component) + (from_com_vec[0] * CoM_component),
					randomNumGen->randNorm(0, random_component) + (from_com_vec[1] * CoM_component),
					randomNumGen->randNorm(0, random_component) + (from_com_vec[2] * CoM_component));
			new_ca_vec += transVec;
		}

		atom_shared_ptr atm = protein->add_new_atom(new_ca_vec,atom_type("CA"),insert_resnum+i);
		//atm->set_b_factor(double(i)+0.1);
		//new_residues[i]->newAtom(CA, new_ca_vec.x, new_ca_vec.y, new_ca_vec.z, 1.0, 0);


	}

	protein->index();

	//chainPtr->renumberExtResNums(1);

	string ext_res_num = c_anchor_res->get_pdb_residue_index();
	trim(ext_res_num);
	const int int_ext_res_num = lexical_cast<int>(ext_res_num);
	int this_res_num = int_ext_res_num-1;
	for (int i = insert_size-1; i >= 0; i--){
		string new_resnum = lexical_cast<string>(this_res_num--);
		(new_residues[i])->set_pdb_residue_index(new_resnum);
	}

	return new_residues;

}



PRODART::POSE::residue_shared_ptr_vector insert_ca_loop_linear_interpolate_minimise(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const double random_component ){

	double CoM_component = 1.0;
	residue_shared_ptr_vector res_vec = insert_ca_loop_linear_interpolate(protein, insert_resnum, rtVec, random_component, CoM_component);
	bool_vector atom_selection(protein->get_all_atom_count(), false);
	for (unsigned int i = 0; i < res_vec.size(); i++){
		const int res_num = res_vec[i]->get_internal_residue_index();
		atom_shared_ptr atm = protein->get_bb_atom(POSE::CA,res_num);
		atom_selection[atm->get_seq_num()] = true;
	}

	ca_minimise_bonds_bumps_angles(protein, atom_selection);
	protein->index();
	atom_selection = bool_vector(protein->get_all_atom_count(), false);
	for (unsigned int i = 0; i < res_vec.size(); i++){
		const int res_num = res_vec[i]->get_internal_residue_index();
		atom_shared_ptr atm = protein->get_bb_atom(POSE::CA,res_num);
		atom_selection[atm->get_seq_num()] = true;
	}
	ca_minimise_bonds_bumps(protein, atom_selection);
	protein->index();
	atom_selection = bool_vector(protein->get_all_atom_count(), false);
	for (unsigned int i = 0; i < res_vec.size(); i++){
		const int res_num = res_vec[i]->get_internal_residue_index();
		atom_shared_ptr atm = protein->get_bb_atom(POSE::CA,res_num);
		atom_selection[atm->get_seq_num()] = true;
	}
	ca_minimise_bonds_bumps_angles(protein, atom_selection);
	/*
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr min_pot = PRODART::POSE::POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("ca_min_default");
	PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(protein);
	minimiser_shared_ptr min_sim =  new_minimiser(min_pot, atom_selection);
	min_sim->make_move(meta_data);
	*/

	return res_vec;
}

PRODART::POSE::residue_shared_ptr_vector insert_ca_loop_linear_interpolate_minimise(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite){

	cout << "start res_count\t" << protein->get_residue_count() << endl;

	for (int i = insert_resnum+num_overwrite-1; i >= insert_resnum; i--){
		protein->delete_residue(i);
	}

	residue_shared_ptr_vector res_vec = insert_ca_loop_linear_interpolate_minimise(protein, insert_resnum, rtVec, 0.01);

	cout  << "end res_count\t" << protein->get_residue_count() << endl;

	return res_vec;
}

PRODART::POSE::residue_shared_ptr_vector remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot){

	// load bb builder
	PRODART::POSE::BB_BUILDER::const_backbone_builder_shared_ptr bb_builder = PRODART::POSE::BB_BUILDER::backbone_builder::Instance();


	residue_shared_ptr_vector res_vec = insert_ca_loop_linear_interpolate_minimise(protein, insert_resnum, rtVec, num_overwrite);

	bool_vector loop_mask(protein->get_residue_count(), false);

	residue_shared_ptr_vector::const_iterator iter;
	for (iter = res_vec.begin(); iter != res_vec.end(); iter++){
		const int resnum = (*iter)->get_internal_residue_index();
		loop_mask[resnum] = true;
	}

	run_ca_loop_anneal(protein,
			loop_mask,
			ca_pot,	//"ca_default",
			high_T_steps,
			anneal_steps,
			low_T_steps,
			start_beta,
			final_beta);

	cout << "remodel_loops: DB point 34235" << endl;

	//lets try new
	CA2MAIN::new_ca2main_minimise(protein, loop_mask, min_steps, bb_min_pot);

	/*
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
	cout << "H_O_pairs:\t" << bb_meta->get_H_O_pair_list().size() << endl;
	cout << "all_bb_pairs:\t" << bb_meta->get_all_bb_pair_list().size() << endl;
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("bb_min_default");
	potentials_energies_map enrg_map = bb_min_pot->get_default_energies_map();
	bb_min_pot->get_energy(bb_meta,enrg_map);
	enrg_map.print_weighted_components(cout);
	cout << endl;

	bb_minimise_geom(protein, loop_mask);
	*/

	return res_vec;
}

PRODART::POSE::residue_shared_ptr_vector remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps,
		const std::string ca_potentials_weight_set,
		const std::string bb_min_potentials_weight_set){
	potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot = pot_factory->make_preset_potentials_container(ca_potentials_weight_set);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = pot_factory->make_preset_potentials_container(bb_min_potentials_weight_set);
	return remodel_loops(protein, insert_resnum, rtVec, num_overwrite, high_T_steps, anneal_steps, low_T_steps, start_beta, final_beta, min_steps, ca_pot , bb_min_pot);
}

bool remodel_loops_multi(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long num_structs,
		std::ostream &output_structs,
		const bool renumber,
		const unsigned long min_steps,
		const std::string ca_potentials_weight_set,
		const std::string bb_min_potentials_weight_set){


	if (num_structs == 0){
		cerr << "remodel_loops_multi: WARNING: 0 structures demanded, returning\n" << endl;
		return true;
	}

	const bool is_const_len = ((unsigned int)num_overwrite == rtVec.size());

	POSE_UTILS::quick_add_HN(protein, false);
	protein->index();

	PRODART::POSE::pose_shared_ptr orig_protein = protein->clone();
	PRODART::POSE::pose_shared_ptr out_protein = protein->clone();


	potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot = pot_factory->make_preset_potentials_container(ca_potentials_weight_set);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = pot_factory->make_preset_potentials_container(bb_min_potentials_weight_set);

	PRODART::POSE::residue_shared_ptr_vector loop_vec = remodel_loops(protein,
			insert_resnum,
			rtVec,
			num_overwrite,
			1000,
			10000,
			1000,
			start_beta,
			final_beta,
			min_steps,
			ca_pot ,
			bb_min_pot);

	if (renumber){
		residue_shared_ptr nres = protein->get_residue(loop_vec[0]->get_chain()->get_first_internal_residue_index());
		string start_ren_str = nres->get_pdb_residue_index();
		trim(start_ren_str);
		const int start_ren = lexical_cast<int>(start_ren_str);
		cout << "renumbering..." << endl;
		protein->renumber_residues(loop_vec[0]->get_chain()->getChainID(), start_ren);

	}

	MISC::inactivate_unset_GLY_CBs(protein);
	cout << "remodel_loops_multi: outputting structure " << 1 << endl;
	//protein->outputPdb(output_structs);

	const int start_resnum = loop_vec[0]->get_internal_residue_index();
	const int end_resnum = loop_vec[loop_vec.size() - 1]->get_internal_residue_index();

	PRODART::POSE::atom_shared_ptr_vector orig_loop;
	PRODART::POSE::atom_shared_ptr_vector moved_loop;


	out_protein = protein->clone();
	if (is_const_len){
		stringstream rmsd_out_string(stringstream::in | stringstream::out);
		orig_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(orig_protein, start_resnum, end_resnum);
		moved_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(protein, start_resnum, end_resnum);
		rmsd_out_string << "loop_RMSD_to_orig:\t" << 1 << "\t" << POSE_UTILS::get_fixed_position_rmsd(orig_loop, moved_loop);
		cout << rmsd_out_string.str() << endl;
		out_protein->add_remark_with_prodart_label("remodel_loops_multi");
		out_protein->add_remark_with_prodart_label(rmsd_out_string.str());
	}
	out_protein->outputPdb(output_structs);
	bool runOK = true;
	unsigned long steps = 0;
	while (steps < (num_structs-1) && runOK){
		runOK = remodel_loops(protein,
				start_resnum,
				end_resnum,
				high_T_steps,
				anneal_steps,
				low_T_steps,
				start_beta,
				final_beta,
				min_steps,
				ca_pot,
				bb_min_pot);
		cout << "remodel_loops_multi: outputting structure " << steps + 2 << endl;

		MISC::inactivate_unset_GLY_CBs(protein);
		out_protein = protein->clone();
		if (is_const_len){
			stringstream rmsd_out_string(stringstream::in | stringstream::out);
			rmsd_out_string << "loop_RMSD_to_orig:\t" << steps + 2 << "\t" << POSE_UTILS::get_fixed_position_rmsd(orig_loop, moved_loop);
			cout << rmsd_out_string.str() << endl;
			out_protein->add_remark_with_prodart_label("remodel_loops_multi");
			out_protein->add_remark_with_prodart_label(rmsd_out_string.str());
		}
		out_protein->outputPdb(output_structs);
		steps++;
	}

	return runOK;
}

//! inserted before insert_resnum
bool remodel_loops_multi_filtered(PRODART::POSE::pose_shared_ptr protein,
		const int insert_resnum,
		const PRODART::POSE::residue_type_vector rtVec,
		const int num_overwrite,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long num_structs,
		std::ostream &output_structs,
		const bool renumber,
		const unsigned long min_steps,
		const std::string ca_potentials_weight_set,
		const std::string bb_min_potentials_weight_set,
		bool no_filter){
	if (num_structs == 0){
		cerr << "remodel_loops_multi: WARNING: 0 structures demanded, returning.\n" << endl;
		return true;
	}

	const bool reset_end_CBs = true;

	MTRand::MTRand_shared_ptr randgen = ENV::get_random_num_gen();


	const bool is_const_len = ((unsigned int)num_overwrite == rtVec.size());

	POSE_UTILS::quick_add_HN(protein, false);
	protein->index();

	PRODART::POSE::pose_shared_ptr orig_protein = protein->clone();
	PRODART::POSE::pose_shared_ptr out_protein = protein->clone();


	potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot = pot_factory->make_preset_potentials_container(ca_potentials_weight_set);
	//PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot_end_cbs = pot_factory->make_preset_potentials_container(ca_potentials_weight_set);
	//POTENTIALS::CA::ca_pseudo_cb_pos_shared_ptr ps_CB_pot = POTENTIALS::CA::new_ca_pseudo_cb_pos();
	//ca_pot_end_cbs->add_potential(ps_CB_pot, 5.0);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot = pot_factory->make_preset_potentials_container(bb_min_potentials_weight_set);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_pure_hs_pot = pot_factory->make_preset_potentials_container("bb_pure_hot_spots");
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_forb_pp_pot = pot_factory->make_preset_potentials_container("bb_forbidden_phi_psi");
	potentials_energies_map hs_enrg_map = bb_pure_hs_pot->get_default_energies_map();

	const int start_resnum = insert_resnum;//loop_vec[0]->get_internal_residue_index();
	const int old_end_resnum = insert_resnum + rtVec.size() - 1 + (rtVec.size() - num_overwrite);
	const int end_resnum = insert_resnum + rtVec.size() - 1;//loop_vec[loop_vec.size() - 1]->get_internal_residue_index();
	if (old_end_resnum != end_resnum){
		cerr << "WARNING: old_end_resnum != end_resnum so may wanna check for consistency" << endl;
	}


	cout << "start_resnum: " << start_resnum << " end_resnum: " << end_resnum << endl;

	/*
	if ( protein->get_bb_atom(POSE::CB, start_resnum-1)->isActiveAndSet()){
		ps_CB_pot->add_rst(start_resnum-1, protein->get_bb_atom(POSE::CB, start_resnum-1)->get_coords());
	}
	if (protein->get_bb_atom(POSE::CB, start_resnum+num_overwrite)->isActiveAndSet()){
		ps_CB_pot->add_rst(end_resnum+1, protein->get_bb_atom(POSE::CB, start_resnum+num_overwrite)->get_coords());
	}
	*/

	PRODART::POSE::residue_shared_ptr_vector loop_vec = remodel_loops(protein,
			insert_resnum,
			rtVec,
			num_overwrite,
			1000,
			10000,
			1000,
			start_beta,
			final_beta,
			min_steps,
			ca_pot ,
			bb_min_pot);


	if (renumber){
		residue_shared_ptr nres = protein->get_residue(loop_vec[0]->get_chain()->get_first_internal_residue_index());
		string start_ren_str = nres->get_pdb_residue_index();
		trim(start_ren_str);
		const int start_ren = lexical_cast<int>(start_ren_str);
		cout << "renumbering..." << endl;
		protein->renumber_residues(loop_vec[0]->get_chain()->getChainID(), start_ren);

	}

	MISC::inactivate_unset_GLY_CBs(protein);



	if (start_resnum != loop_vec[0]->get_internal_residue_index()){
		cerr << "EROR in start_resnum\t" <<  start_resnum << "\t" << loop_vec[0]->get_internal_residue_index() << endl;
	}
	if (end_resnum != loop_vec[loop_vec.size() - 1]->get_internal_residue_index()){
		cerr << "EROR in end_resnum\t" <<  end_resnum << "\t" << loop_vec[loop_vec.size() - 1]->get_internal_residue_index() << endl;
	}
	bool_vector loop_mask(protein->get_residue_count(), false);
	residue_shared_ptr_vector::const_iterator iter;
	for (int i = start_resnum; i <= end_resnum; i++){
		loop_mask[i] = true;
	}

	int start_resnum_m1 = start_resnum;
	int end_resnum_p1 = end_resnum;

	bool_vector loop_mask_ext1 = loop_mask;
	if (start_resnum-1 > 0){
		loop_mask_ext1[start_resnum-1] = true;
		start_resnum_m1 = start_resnum - 1;
	}
	if (end_resnum+1 < int(loop_mask.size())){
		loop_mask_ext1[end_resnum+1] = true;
		end_resnum_p1 = end_resnum + 1;
	}

	POSE::atom_shared_ptr_vector ca_loop_atms;
	for (int i = start_resnum_m1; i <= end_resnum_p1; i++){
		ca_loop_atms.push_back(protein->get_bb_atom(POSE::CA, i));
	}

	const atom_shared_ptr_bool_map orig_act_map = get_activity_map(protein);
	inactive_all_outside_radius(protein, ca_loop_atms, 10.0);
	const atom_shared_ptr_bool_map masked_act_map = get_activity_map(protein);

	bool_vector hot_spot_vec(protein->get_residue_count(), false);

	const double max_forbidden_pp = std::max(0.1 * double(end_resnum - start_resnum + 1), 0.0);
	const double max_forbidden_f3 = std::max(0.7 * double(end_resnum - start_resnum + 1), 2.0);
	const double max_CB_dev = 0.7;

	cout << "max_forbidden_pp " << max_forbidden_pp << endl;
	//cout << "max_forbidden_f3 " << max_forbidden_f3 << endl;

	PRODART::POSE::atom_shared_ptr_vector orig_loop;
	PRODART::POSE::atom_shared_ptr_vector moved_loop;
	const PRODART::POSE::atom_shared_ptr_vector atom_selection = CA2MAIN::get_loop_moved_atoms(protein, start_resnum, end_resnum);
	LOOP::atom_shared_ptr_vector3d_map backup_coords = LOOP::get_coords( atom_selection);

	out_protein = protein->clone();
	if (is_const_len){
		stringstream rmsd_out_string(stringstream::in | stringstream::out);
		orig_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(orig_protein, start_resnum, end_resnum);
		moved_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(protein, start_resnum, end_resnum);
		rmsd_out_string << "loop_RMSD_to_orig:\t" << 0 << "\t" << POSE_UTILS::get_fixed_position_rmsd(orig_loop, moved_loop);
		cout << rmsd_out_string.str() << endl;
		out_protein->add_remark_with_prodart_label("remodel_loops_multi");
		out_protein->add_remark_with_prodart_label(rmsd_out_string.str());
	}
	//cout << "remodel_loops_multi: outputting structure " << 0 << endl;
	//out_protein->outputPdb(output_structs);

	bb_pose_meta_shared_ptr hs_pose_meta = new_bb_pose_meta(protein);
	//double hs_energy =
	bb_pure_hs_pot->get_energy_residue_loop(hs_pose_meta,hs_enrg_map, loop_mask_ext1);
	double acc_forb = hs_enrg_map[potentials_name("bb_forbidden_phi_psi")].energy;
	double acc_f3 = hs_enrg_map[potentials_name("bb_frag3_mq")].energy;
	bb_pure_hs_pot->get_residue_hot_spots(hs_pose_meta, hot_spot_vec);
	double hs_count = 0;
	for (unsigned int it = 0; it != hot_spot_vec.size(); it++){
		if (hot_spot_vec[it] && loop_mask[it] ) hs_count += 1;
	}

	//bool_vector hs_loop_mask = loop_mask;

	double prev_sum_score = 99999999.9;
	unsigned long tries = 0;
	bool runOK = true;
	unsigned long steps = 0;
	while (steps < (num_structs) && runOK){

		hs_pose_meta = new_bb_pose_meta(protein);
		bb_forb_pp_pot->get_residue_hot_spots(hs_pose_meta, hot_spot_vec);
		double this_hs_count = 0;
		for (unsigned int it = 1; it != hot_spot_vec.size()-1; it++){
			if (loop_mask[it] && (hot_spot_vec[it] || hot_spot_vec[it-1] || hot_spot_vec[it+1]) ) this_hs_count += 1;
		}
		bool_vector hs_loop_mask = loop_mask;
		cout << "hscount: " << this_hs_count << endl;;
		if (this_hs_count > 0.0){
			//cout << "db1" << endl;
			const double p_hs_move = randgen->rand(1.0);
			if (p_hs_move < 0.2){
				//cout << "db2" << endl;
				cout <<  "0";
				for (unsigned int hi = 1; hi < hot_spot_vec.size() -1; hi++){
					hs_loop_mask[hi] = loop_mask[hi] && (hot_spot_vec[hi] || hot_spot_vec[hi-1] || hot_spot_vec[hi+1]);
					cout << hs_loop_mask[hi];
				}
				cout <<  "0" << endl;
			}
		}

		if ((tries >= 500) && (tries % 50 == 0)){
			//LOOP::ca_minimise_bonds_bumps(protein, );
			cout << "minimising CAs again" << endl;
			bool_vector ca_atom_selection(protein->get_all_atom_count(), false);
			for (unsigned int i = 0; i < loop_mask.size(); i++){
				if (loop_mask[i]){
					const int res_num = i;
					atom_shared_ptr atm = protein->get_bb_atom(POSE::CA,res_num);
					ca_atom_selection[atm->get_seq_num()] = true;
					if (tries % 100 == 0){
						vector3d vec = atm->get_coords();
						//cout << "randomising CAs again" << endl;
						vec += vector3d(randgen->rand(0.25), randgen->rand(0.25), randgen->rand(0.25));
					}
				}
			}
			ca_minimise_bonds_bumps(protein, ca_atom_selection);
		}

		const bool long_reset = ((steps % 1000 == 0) && (steps !=0) && (tries % 100 == 0)  ) ? true : false;

		if (!long_reset){
		run_ca_loop_anneal(protein,
				hs_loop_mask, //loop_mask,
				ca_pot,	//"ca_default",
				high_T_steps,
				anneal_steps,
				low_T_steps,
				start_beta,
				final_beta);
		}
		else {
			cout << "doing long_reset\n";
			run_ca_loop_anneal(protein,
					hs_loop_mask, //loop_mask,
					ca_pot,	//"ca_default",
					1000,
					10000,
					1000,
					start_beta,
					final_beta);
		}

		//lets try new method
		cout << "DB point" << endl;
		CA2MAIN::new_ca2main(protein, loop_mask);
		MISC::inactivate_unset_GLY_CBs(protein);

		if (reset_end_CBs == true){
			if (protein->get_bb_atom(CB, start_resnum_m1)->isActiveAndSet()){
				protein->get_bb_atom(CB, start_resnum_m1)->set_coords(POSE_UTILS::get_rough_ideal_CB(protein, start_resnum_m1));
			}
			if (protein->get_bb_atom(CB, end_resnum_p1)->isActiveAndSet()){
				protein->get_bb_atom(CB, end_resnum_p1)->set_coords(POSE_UTILS::get_rough_ideal_CB(protein, end_resnum_p1));
			}
		}

		double start_cb_from_ideal = protein->get_bb_atom(CB, start_resnum_m1)->isActiveAndSet()
				? (protein->get_bb_atom(CB, start_resnum_m1)->get_coords() - POSE_UTILS::get_rough_ideal_CB(protein, start_resnum_m1)).mod()
				: 0.0;

		double end_cb_from_ideal = protein->get_bb_atom(CB, end_resnum_p1)->isActiveAndSet()
				? (protein->get_bb_atom(CB, end_resnum_p1)->get_coords() - POSE_UTILS::get_rough_ideal_CB(protein, end_resnum_p1)).mod()
				: 0.0;

		cout << "CBs dist from ideal:\t" << start_cb_from_ideal << "\t" << end_cb_from_ideal << endl;

		hs_pose_meta = new_bb_pose_meta(protein);
		//double this_hs_energy =
		bb_pure_hs_pot->get_energy_residue_loop(hs_pose_meta,hs_enrg_map, loop_mask_ext1);
		bb_pure_hs_pot->get_residue_hot_spots(hs_pose_meta, hot_spot_vec);
		this_hs_count = 0;
		for (unsigned int it = 0; it != hot_spot_vec.size(); it++){
			if (hot_spot_vec[it] && loop_mask[it] ) this_hs_count += 1;
		}
		const double forb = hs_enrg_map[potentials_name("bb_forbidden_phi_psi")].energy;
		const double f3 = hs_enrg_map[potentials_name("bb_frag3_mq")].energy;

		const double this_sum_score = (5.0*forb) +  f3 + (5.0 * start_cb_from_ideal) + (5.0 * end_cb_from_ideal);

		cout << "sum_score:\t" << this_sum_score << endl;

		if ( no_filter ||	(( (start_cb_from_ideal < max_CB_dev) && (end_cb_from_ideal < max_CB_dev))
				&& (((f3 <= std::max(acc_f3, max_forbidden_f3)) && (forb <= std::max(acc_forb, max_forbidden_pp) ))
				|| (forb <= max_forbidden_pp && f3 <= max_forbidden_f3)) ) ){

			CA2MAIN::bb_minimise(protein, loop_mask, bb_min_pot, min_steps);

			cout << "remodel_loops_multi: outputting structure " << steps + 1 << endl;

			//MISC::inactivate_unset_GLY_CBs(protein);
			restore_activity_map(orig_act_map);
			out_protein = protein->clone();
			if (is_const_len){
				stringstream rmsd_out_string(stringstream::in | stringstream::out);
				rmsd_out_string << "loop_RMSD_to_orig:\t" << steps + 1 << "\t" << POSE_UTILS::get_fixed_position_rmsd(orig_loop, moved_loop);
				cout << rmsd_out_string.str() << endl;
				out_protein->add_remark_with_prodart_label("remodel_loops_multi");
				out_protein->add_remark_with_prodart_label(rmsd_out_string.str());
			}
			hs_enrg_map.print_headers(cout);
			hs_enrg_map.print_unweighted_components(cout);
			out_protein->outputPdb(output_structs);
			restore_activity_map(masked_act_map);
			//hs_energy = this_hs_energy;
			hs_count = this_hs_count;
			acc_forb = forb;
			acc_f3 = f3;
			backup_coords = get_coords(atom_selection);
			tries = 0;
			steps++;
		}
		else {
			const double prob = randgen->rand(1.0);

			if ((this_sum_score <= prev_sum_score)){
				backup_coords = get_coords(atom_selection);
			}
			else if (long_reset){
				backup_coords = get_coords(atom_selection);
			}
			else if ((prob < 0.6) ){
				restore_coords(backup_coords);
			}
			else {
				backup_coords = get_coords(atom_selection);
			}
			hs_enrg_map.print_headers(cout);
			hs_enrg_map.print_unweighted_components(cout);
			cout << "Rejected" << endl;
			tries++;
		}

		prev_sum_score = this_sum_score;

	}

	return runOK;
}

double get_fiser_bb_loop_rmsd(PRODART::POSE::pose_shared_ptr ref,
		PRODART::POSE::pose_shared_ptr protein,
		const int start_res,
		const int length){

	const int end_resnum = start_res + length - 1;

	PRODART::POSE::atom_shared_ptr_vector orig_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(ref, start_res, end_resnum);
	PRODART::POSE::atom_shared_ptr_vector moved_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(protein, start_res, end_resnum);
	return POSE_UTILS::get_fixed_position_rmsd(orig_loop, moved_loop);


}

double get_fiser_bb_loop_rmsd_corrected(PRODART::POSE::pose_shared_ptr ref,
		PRODART::POSE::pose_shared_ptr protein,
		const int start_res,
		const int length){

	const int end_resnum = start_res + length - 1;

	PRODART::POSE::atom_shared_ptr_vector orig_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(ref, start_res+1, end_resnum-1);
	PRODART::POSE::atom_shared_ptr_vector moved_loop = CA2MAIN::get_loop_moved_N_CA_C_O_atoms(protein, start_res+1, end_resnum-1);
	return POSE_UTILS::get_fixed_position_rmsd(orig_loop, moved_loop);


}

bool run_ca_loop_anneal(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta){
	//potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	//PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot = pot_factory->make_preset_potentials_container(potentials_weight_set);


	PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(protein);

	// TODO replace this with mover factory
	/*
	PRODART::POSE::MOVERS::move_set_shared_ptr dih_move_set = PRODART::POSE::MOVERS::CA::ca_single_dihedral_uni_dist_move_set_factory(meta_data, 0.0038);
	PRODART::POSE::MOVERS::move_set_shared_ptr ang_move_set = PRODART::POSE::MOVERS::CA::ca_single_angle_uni_dist_move_set_factory(meta_data, 0.0038);
	PRODART::POSE::MOVERS::move_set_shared_ptr pinch_move_set = PRODART::POSE::MOVERS::CA::ca_local_angle_pinch_uni_dist_move_set_factory(meta_data, 0.0038);
	PRODART::POSE::MOVERS::move_set_shared_ptr bond_move_set = PRODART::POSE::MOVERS::CA::ca_local_bond_uni_dist_move_set_factory(meta_data, 0.01);
	PRODART::POSE::MOVERS::move_set_shared_ptr crank_move_set = PRODART::POSE::MOVERS::CA::ca_single_crankshaft_uni_dist_move_set_factory(meta_data,
			0.01,
			2,
			8);
	PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = PRODART::POSE::MOVERS::new_move_set();
	overall_move_set->add_move(dih_move_set, 1.0);
	overall_move_set->add_move(ang_move_set, 1.0);
	overall_move_set->add_move(bond_move_set, 1.0);
	overall_move_set->add_move(crank_move_set, 1.0);
	overall_move_set->add_move(pinch_move_set, 1.0);
	overall_move_set->propagate_rand_num_gen(PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id()));
	*/

	PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = MOVERS::move_set_factory::Instance()->make_preset_loop_move_set_by_residue(meta_data,
			"ca_standard",
			PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id()),
			loop_mask);


	PRODART::POSE::SIM::monte_carlo_shared_ptr mc_sim = PRODART::POSE::SIM::new_monte_carlo(pot,
			overall_move_set,
			PRODART::POSE::SIM::CA::new_ca_motif_anneal_protocol(high_T_steps, //new_simple_ca_sim_anneal_protocol
					anneal_steps,
					low_T_steps,
					start_beta,
					final_beta,
					PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id())));


	cout << "running anneal" << endl;

	mc_sim->make_move(meta_data);

	return true;
}

bool run_ca_loop_anneal(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const std::string potentials_weight_set,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta){

	potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot = pot_factory->make_preset_potentials_container(potentials_weight_set);


	PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(protein);

	// TODO replace this with mover factory
	/*
	PRODART::POSE::MOVERS::move_set_shared_ptr dih_move_set = PRODART::POSE::MOVERS::CA::ca_single_dihedral_uni_dist_move_set_factory(meta_data, 0.0038);
	PRODART::POSE::MOVERS::move_set_shared_ptr ang_move_set = PRODART::POSE::MOVERS::CA::ca_single_angle_uni_dist_move_set_factory(meta_data, 0.0038);
	PRODART::POSE::MOVERS::move_set_shared_ptr pinch_move_set = PRODART::POSE::MOVERS::CA::ca_local_angle_pinch_uni_dist_move_set_factory(meta_data, 0.0038);
	PRODART::POSE::MOVERS::move_set_shared_ptr bond_move_set = PRODART::POSE::MOVERS::CA::ca_local_bond_uni_dist_move_set_factory(meta_data, 0.01);
	PRODART::POSE::MOVERS::move_set_shared_ptr crank_move_set = PRODART::POSE::MOVERS::CA::ca_single_crankshaft_uni_dist_move_set_factory(meta_data,
			0.01,
			2,
			8);
	PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = PRODART::POSE::MOVERS::new_move_set();
	overall_move_set->add_move(dih_move_set, 1.0);
	overall_move_set->add_move(ang_move_set, 1.0);
	overall_move_set->add_move(bond_move_set, 1.0);
	overall_move_set->add_move(crank_move_set, 1.0);
	overall_move_set->add_move(pinch_move_set, 1.0);
	overall_move_set->propagate_rand_num_gen(PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id()));
	*/

	PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = MOVERS::move_set_factory::Instance()->make_preset_loop_move_set_by_residue(meta_data,
			"ca_standard",
			PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id()),
			loop_mask);


	PRODART::POSE::SIM::monte_carlo_shared_ptr mc_sim = PRODART::POSE::SIM::new_monte_carlo(pot,
			overall_move_set,
			PRODART::POSE::SIM::CA::new_ca_motif_anneal_protocol(high_T_steps, //new_simple_ca_sim_anneal_protocol
					anneal_steps,
					low_T_steps,
					start_beta,
					final_beta,
					PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id())));


	cout << "running anneal" << endl;

	mc_sim->make_move(meta_data);

	return true;
}

bool ca_minimise_bonds_bumps(PRODART::POSE::pose_shared_ptr protein){
	protein->index();
	bool_vector atom_selection(protein->get_all_atom_count(), false);
	for (int i = 0; i < protein->get_residue_count(); i++){
		atom_shared_ptr at = protein->get_bb_atom(POSE::CA, i);
		if (at){
			if (at->isActiveAndSet()){
				atom_selection[at->get_seq_num()] = true;
			}
		}

	}
	return ca_minimise_bonds_bumps(protein, atom_selection);
}

bool ca_minimise_bonds_bumps(PRODART::POSE::pose_shared_ptr protein,
		bool_vector& atom_selection){
	PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(protein);

	const int resCount = protein->get_residue_count();

	if (resCount > 1){

		//residue_shared_ptr prev_res;

		for (int i = 0; i < resCount-1; i++){
			residue_shared_ptr this_res = protein->get_residue(i);
			residue_shared_ptr this_res_p1 = protein->get_residue(i+1);
			if (this_res_p1 && this_res){

				if (this_res_p1->get_chain() == this_res->get_chain()){

					atom_shared_ptr at1 = protein->get_bb_atom(POSE::CA, i);
					atom_shared_ptr at2 = protein->get_bb_atom(POSE::CA, i+1);
					if (at1 && at2){
						// add bond

						const int at1_num = at1->get_seq_num();
						const int at2_num = at2->get_seq_num();

						if (atom_selection[at1_num] || atom_selection[at2_num]){
							//cout << "added bond " << at1_num << " " << at2_num << endl;
							meta_data->add_dist_harmonic_restraint(at1->get_seq_num(),
									at2->get_seq_num(),
									3.80581,
									890.333);
						}

					}


				}

			}
			//prev_res = this_res;
		}
	}


	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr min_pot = PRODART::POSE::POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("ca_min_default");
	potentials_energies_map enrg_map = min_pot->get_default_energies_map();
	min_pot->get_energy(meta_data,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	minimiser_interface_shared_ptr min_sim =  new_minimiser(min_pot, atom_selection);
	bool result = min_sim->make_move(meta_data);
	min_pot->get_energy(meta_data,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	return result;
}

bool ca_minimise_bonds_bumps_angles(PRODART::POSE::pose_shared_ptr protein,
		bool_vector& atom_selection,
		bool auto_load_restraints){
	PRODART::POSE::META::ca_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_ca_pose_meta(protein,auto_load_restraints);

	const int resCount = protein->get_residue_count();

	if (resCount > 1){

		//residue_shared_ptr prev_res;

		for (int i = 0; i < resCount-1; i++){
			residue_shared_ptr this_res = protein->get_residue(i);
			residue_shared_ptr this_res_p1 = protein->get_residue(i+1);
			if (this_res_p1 && this_res){

				if (this_res_p1->get_chain() == this_res->get_chain()){

					atom_shared_ptr at1 = protein->get_bb_atom(POSE::CA, i);
					atom_shared_ptr at2 = protein->get_bb_atom(POSE::CA, i+1);
					if (at1 && at2){
						// add bond

						const int at1_num = at1->get_seq_num();
						const int at2_num = at2->get_seq_num();

						if (atom_selection[at1_num] || atom_selection[at2_num]){
							meta_data->add_dist_harmonic_restraint(at1->get_seq_num(),
									at2->get_seq_num(),
									3.80581,
									890.333);
						}

					}


				}

			}
			//prev_res = this_res;
		}

		for (int i = 0; i < resCount-2; i++){
			residue_shared_ptr this_res = protein->get_residue(i);
			residue_shared_ptr this_res_p1 = protein->get_residue(i+1);
			residue_shared_ptr this_res_p2 = protein->get_residue(i+2);
			if (this_res_p1 && this_res && this_res_p2 ){

				if ((this_res_p1->get_chain() == this_res->get_chain())
						&& (this_res_p1->get_chain() == this_res_p2->get_chain())){

					atom_shared_ptr at1 = protein->get_bb_atom(POSE::CA, i);
					atom_shared_ptr at2 = protein->get_bb_atom(POSE::CA, i+1);
					atom_shared_ptr at3 = protein->get_bb_atom(POSE::CA, i+2);
					if (at1 && at2 && at3){
						// add bond

						const int at1_num = at1->get_seq_num();
						const int at2_num = at2->get_seq_num();
						const int at3_num = at3->get_seq_num();

						if (atom_selection[at1_num] || atom_selection[at2_num] || atom_selection[at3_num] ){
							meta_data->add_angle_harmonic_restraint(at1->get_seq_num(),
									at2->get_seq_num(),
									at3->get_seq_num(),
									UTILS::degrees_to_radians(110),
									50.0);
						}

					}


				}

			}
			//prev_res = this_res;
		}
	}


	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr min_pot = PRODART::POSE::POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("ca_min_default");
	potentials_energies_map enrg_map = min_pot->get_default_energies_map();
	min_pot->get_energy(meta_data,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	minimiser_interface_shared_ptr min_sim =  new_minimiser(min_pot, atom_selection);
	bool result = min_sim->make_move(meta_data);
	min_pot->get_energy(meta_data,enrg_map);
	enrg_map.print_headers(cout);
	enrg_map.print_weighted_components(cout);
	return result;
}



bool remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ca_pot,
		PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_min_pot){
	// load bb builder
	PRODART::POSE::BB_BUILDER::const_backbone_builder_shared_ptr bb_builder = PRODART::POSE::BB_BUILDER::backbone_builder::Instance();

	bool_vector loop_mask(protein->get_residue_count(), false);

	residue_shared_ptr_vector::const_iterator iter;
	for (int i = start_resnum; i <= end_resnum; i++){
		loop_mask[i] = true;
	}

	run_ca_loop_anneal(protein,
			loop_mask,
			ca_pot,
			high_T_steps,
			anneal_steps,
			low_T_steps,
			start_beta,
			final_beta);

	CA2MAIN::simple_ca2main_minimise(protein, loop_mask, min_steps, bb_min_pot);



	return true;
}

bool remodel_loops(PRODART::POSE::pose_shared_ptr protein,
		const int start_resnum,
		const int end_resnum,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const unsigned long min_steps,
		const std::string ca_potentials_weight_set,
		const std::string bb_min_potentials_weight_set){

	// load bb builder
	PRODART::POSE::BB_BUILDER::const_backbone_builder_shared_ptr bb_builder = PRODART::POSE::BB_BUILDER::backbone_builder::Instance();

	bool_vector loop_mask(protein->get_residue_count(), false);

	residue_shared_ptr_vector::const_iterator iter;
	for (int i = start_resnum; i <= end_resnum; i++){
		loop_mask[i] = true;
	}

	run_ca_loop_anneal(protein,
			loop_mask,
			ca_potentials_weight_set,
			high_T_steps,
			anneal_steps,
			low_T_steps,
			start_beta,
			final_beta);

	CA2MAIN::simple_ca2main_minimise(protein, loop_mask, min_steps, bb_min_potentials_weight_set);



	return true;

}

bool_vector get_bb_hot_spots(PRODART::POSE::pose_shared_ptr protein){
	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_hot_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("bb_hot_spots_default");
	bool_vector hot_spots(protein->get_residue_count(), false);
	bb_hot_pot->get_residue_hot_spots(bb_meta, hot_spots);

	return hot_spots;
}

void print_bb_hot_spots(PRODART::POSE::pose_shared_ptr protein, std::ostream& output){
	bool_vector hot_spots = get_bb_hot_spots(protein);
	for (unsigned int i = 0; i < hot_spots.size(); i++){
		if(hot_spots[i]){
			const_residue_shared_ptr res = protein->get_residue(i);
			output << res->get_pdb_residue_index() << "\t"
					<< res->get_chain()->getChainID() << "\t"
					<< res->get_type().get_label() << "\t"
					<< endl;
		}
	}



	bb_pose_meta_shared_ptr bb_meta = POSE::META::new_bb_pose_meta(protein);
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr bb_hot_pot = POTENTIALS::potentials_factory::Instance()->make_preset_potentials_container("bb_hot_spots_default");
	potentials_energies_map energies = bb_hot_pot->get_default_energies_map();
	bb_hot_pot->get_energy(bb_meta, energies);
	energies.print_headers(output);
	energies.print_weighted_components(output);

	bool_vector::const_iterator iter;
	for (iter = hot_spots.begin(); iter != hot_spots.end(); iter++){
		output << *iter;
	}
	output << endl;
	cout << bb_meta->get_sec_struct() << endl;

}
bool fixed_bb_motif_anneal(PRODART::POSE::pose_shared_ptr protein,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_wt,
		const double final_wt){
	POSE_UTILS::quick_add_HN(protein, false);
	MISC::inactivate_unset_GLY_CBs(protein);
	protein->index();

	const bool_vector as_site_mask(protein->get_residue_count(), true);
	const bool_vector loop_mask(protein->get_residue_count(), true);

	potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	POTENTIALS::as_geom_pot_shared_ptr as_pot = boost::static_pointer_cast<POTENTIALS::as_geom_pot, POTENTIALS::potential_interface>(PRODART::POSE::POTENTIALS::potentials_factory::Instance()->get_potential(POTENTIALS::potentials_name("as_geom_pot")));
	as_pot->set_loop_mask(as_site_mask);
	as_pot->random_assign_motif(protein);


	//as_pot->do_exaustive_bb_motif_search(protein);

	potentials_name_double_map wts;
	wts[potentials_name("as_geom_pot")] = 1.0;
	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot = pot_factory->make_potentials_container(wts);
	PRODART::POSE::META::bb_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_bb_pose_meta(protein);

	PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = MOVERS::move_set_factory::Instance()->make_preset_loop_move_set_by_residue(meta_data,
			"bb_composite_ca",
			PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id()),
			loop_mask);

	overall_move_set->propagate_bb_potential_preset("bb_min_default_as_geom_pot");
	overall_move_set->propagate_ca_potential_preset("ca_default_as_geom_pot");


	PRODART::POSE::SIM::monte_carlo_shared_ptr mc_sim = PRODART::POSE::SIM::new_monte_carlo(pot,
			overall_move_set,
			PRODART::POSE::SIM::BB::new_fixed_bb_motif_anneal(high_T_steps, anneal_steps, low_T_steps, start_wt, final_wt,
			PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id())));

	cout << "running mc" << endl;


	const bool result =  mc_sim->make_move(meta_data);

	const double rmsd = as_pot->get_bb_rmsd(protein);

	cout << "transfering sidechains:" << endl;
	as_pot->transfer_sidechains(meta_data);
	cout << "transfering non-peptides:" << endl;
	as_pot->transfer_non_peptides(meta_data);
	cout << "done." << endl;

	string remark("as_geom_pot: bb_RMSD:\t");
	remark.append(lexical_cast<string>(rmsd));
	protein->add_remark_with_prodart_label(remark);

	mc_sim->get_state().print_summary(cout);

	cout << protein->get_chain_count() << endl;



	return result;
}

bool fixed_bb_motif_brute_force_search(PRODART::POSE::pose_shared_ptr protein){
	//potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	POTENTIALS::as_geom_pot_shared_ptr as_pot = boost::static_pointer_cast<POTENTIALS::as_geom_pot, POTENTIALS::potential_interface>(PRODART::POSE::POTENTIALS::potentials_factory::Instance()->get_potential(POTENTIALS::potentials_name("as_geom_pot")));
	const bool_vector as_site_mask(protein->get_residue_count(), true);
	const bool_vector loop_mask(protein->get_residue_count(), true);
	as_pot->set_loop_mask(as_site_mask);
	as_pot->random_assign_motif(protein);
	PRODART::POSE::META::bb_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_bb_pose_meta(protein);

	const double best_rmsd = as_pot->do_exaustive_bb_motif_search(protein);

	cout << "BEST_RMSD: " << best_rmsd << endl;

	//const double rmsd = as_pot->get_bb_rmsd(protein);

	as_pot->print_assign_info(cout);

	cout << "transfering sidechains:" << endl;
	as_pot->transfer_sidechains(meta_data);
	cout << "transfering non-peptides:" << endl;
	as_pot->transfer_non_peptides(meta_data);
	cout << "done." << endl;

	string remark("as_geom_pot: bb_RMSD:\t");
	remark.append(lexical_cast<string>(best_rmsd));
	protein->add_remark_with_prodart_label(remark);

	cout << protein->get_chain_count() << endl;

	return true;
}

bool ca_bb_multi_level_loop_design(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const bool_vector& as_site_mask,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta){


	POSE_UTILS::quick_add_HN(protein, false);
	MISC::inactivate_unset_GLY_CBs(protein);
	protein->index();
	//CA2MAIN::simple_ca2main_minimise(in_pdb);

	if (loop_mask.size() != (unsigned int)protein->get_residue_count()){
		std::cerr << "ca_bb_multi_level_loop_design: ERROR: loop_mask different size to protein" << endl;
		return false;
	}

	potentials_factory *pot_factory = PRODART::POSE::POTENTIALS::potentials_factory::Instance();
	POTENTIALS::as_geom_pot_shared_ptr as_pot = boost::static_pointer_cast<POTENTIALS::as_geom_pot, POTENTIALS::potential_interface>(PRODART::POSE::POTENTIALS::potentials_factory::Instance()->get_potential(POTENTIALS::potentials_name("as_geom_pot")));
	as_pot->set_loop_mask(as_site_mask);
	as_pot->random_assign_motif(protein);

	PRODART::POSE::POTENTIALS::potentials_container_shared_ptr pot = pot_factory->make_preset_potentials_container("bb_hot_spots_default_as_geom_pot");
	PRODART::POSE::META::bb_pose_meta_shared_ptr meta_data = PRODART::POSE::META::new_bb_pose_meta(protein);

	PRODART::POSE::MOVERS::move_set_shared_ptr overall_move_set = MOVERS::move_set_factory::Instance()->make_preset_loop_move_set_by_residue(meta_data,
			"bb_composite_ca",
			PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id()),
			loop_mask);

	overall_move_set->propagate_bb_potential_preset("bb_min_default_as_geom_pot");
	overall_move_set->propagate_ca_potential_preset("ca_default_as_geom_pot");


	PRODART::POSE::SIM::monte_carlo_shared_ptr mc_sim = PRODART::POSE::SIM::new_monte_carlo(pot,
			overall_move_set,
			PRODART::POSE::SIM::BB::new_bb_motif_anneal_protocol(high_T_steps, anneal_steps, low_T_steps, start_beta, final_beta,
			PRODART::ENV::prodart_env::Instance()->get_random_num_gen(boost::this_thread::get_id())));

	cout << "running mc" << endl;


	const bool result =  mc_sim->make_move(meta_data);

	const double rmsd = as_pot->get_bb_rmsd(protein);

	cout << "transfering sidechains:" << endl;
	as_pot->transfer_sidechains(meta_data);
	cout << "transfering non-peptides:" << endl;
	as_pot->transfer_non_peptides(meta_data);
	cout << "done." << endl;

	string remark("as_geom_pot: bb_RMSD:\t");
	remark.append(lexical_cast<string>(rmsd));
	protein->add_remark_with_prodart_label(remark);

	cout << protein->get_chain_count() << endl;

	return result;
}

bool ca_bb_multi_level_loop_design(PRODART::POSE::pose_shared_ptr protein,
		const bool_vector& loop_mask,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta){
	return ca_bb_multi_level_loop_design( protein,
			loop_mask,
			loop_mask,
			high_T_steps,
			anneal_steps,
			low_T_steps,
			start_beta,
			final_beta);
}

atom_shared_ptr_bool_map get_activity_map(PRODART::POSE::pose_shared_ptr protein){

	const int num_ats = protein->get_all_atom_count();
	atom_shared_ptr_bool_map rtn_map;

	for (int i = 0 ; i < num_ats; i++){
		PRODART::POSE::atom_shared_ptr atm = protein->get_atom(i);
		if (atm->isSet()){
			rtn_map[atm] = atm->isActive();
		}
	}

	return rtn_map;

}

void restore_activity_map(const atom_shared_ptr_bool_map& vec){

	atom_shared_ptr_bool_map::const_iterator it;

	for (it = vec.begin(); it != vec.end(); it++){
		it->first->setActive(it->second);
	}


}


void inactive_all_outside_radius(PRODART::POSE::pose_shared_ptr protein, const POSE::atom_shared_ptr_vector loop, const double radius ){


	const int num_ats = protein->get_all_atom_count();
	//const int num_loop_ats = loop.size();
	atom_shared_ptr_bool_map rtn_map;

	for (int i = 0 ; i < num_ats; i++){
		PRODART::POSE::atom_shared_ptr atm = protein->get_atom(i);
		const vector3d at_vec = atm->get_coords();
		bool is_within = false;
		if (atm->isSet()){
			for (atom_shared_ptr_vector::const_iterator it = loop.begin(); it != loop.end(); it++){
				const vector3d loop_vec = (*it)->get_coords();
				if ((at_vec - loop_vec).mod() <= radius){
					is_within = true;
				}
			}
		}
		if (is_within == false){
			atm->setActive(false);
		}
	}

}


}
}
}







