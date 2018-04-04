//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose_basic_kine.cpp
 *
 *  Created on: 23 Feb 2010
 *      Author: jmacdona
 */
#include "pose_basic_kine.h"


using namespace PRODART::POSE;
using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::pow;
using std::sqrt;
using std::asin;
using std::string;

using boost::lexical_cast;
using boost::split;
using boost::is_any_of;
using boost::trim;

using std::ios;
using std::setiosflags;
using std::resetiosflags;
using std::setprecision;
using std::setw;

using PRODART::UTILS::vector3d;
using namespace PRODART::UTILS;
using namespace PRODART::POSE;


namespace PRODART {
namespace POSE_UTILS {
namespace KINE {





void translate_fa(const PRODART::POSE::pose_shared_ptr protein,
		const PRODART::UTILS::vector3d& trans){
	const int bb_atom_count = protein->get_bb_atom_count();

	for (int i = 0; i < bb_atom_count; i++){
		vector3d vec = protein->get_bb_atom_coords(i);
		vec = vec + trans;
		protein->set_bb_atom_coords(vec, i);
	}

	const int res_count = protein->get_residue_count();

	for (int i = 0; i < res_count; i++){
		residue_shared_ptr currRes = protein->get_residue(i);
		sidechain_shared_ptr currSC = currRes->get_sidechain();
		const int sc_atom_count = currSC->get_atom_count();
		for (int j = 0; j < sc_atom_count; j++){
			vector3d vec = currSC->get_atom(j)->get_coords();
			vec = vec + trans;
			currSC->get_atom(j)->set_coords(vec);
		}
	}
}

void translate_fa(const PRODART::POSE::atom_shared_ptr_vector atoms_to_move,
		const PRODART::UTILS::vector3d& trans){
	const int atom_count = atoms_to_move.size();
	for (int i = 0; i < atom_count; i++){
		atom_shared_ptr ptr = atoms_to_move[i];
		vector3d vec = ptr->get_coords();
		vec = vec + trans;
		ptr->set_coords(vec);
	}
}

void apply_rotation_matrix_fa(const PRODART::POSE::pose_shared_ptr protein,
		const PRODART::UTILS::rot_matrix& rot){

	const int bb_atom_count = protein->get_bb_atom_count();

	for (int i = 0; i < bb_atom_count; i++){
		vector3d vec = protein->get_bb_atom_coords(i);
		PRODART::UTILS::matrix_rotate(vec, rot);
		protein->set_bb_atom_coords(vec, i);
	}

	const int res_count = protein->get_residue_count();

	for (int i = 0; i < res_count; i++){
		residue_shared_ptr currRes = protein->get_residue(i);
		sidechain_shared_ptr currSC = currRes->get_sidechain();
		const int sc_atom_count = currSC->get_atom_count();
		for (int j = 0; j < sc_atom_count; j++){
			vector3d vec = currSC->get_atom(j)->get_coords();
			PRODART::UTILS::matrix_rotate(vec, rot);
			currSC->get_atom(j)->set_coords(vec);
		}
	}

}

void apply_rotation_matrix_fa(const PRODART::POSE::atom_shared_ptr_vector atoms_to_move,
		const PRODART::UTILS::rot_matrix& rot){
	const int atom_count = atoms_to_move.size();
	for (int i = 0; i < atom_count; i++){
		atom_shared_ptr ptr = atoms_to_move[i];
		vector3d vec = ptr->get_coords();
		PRODART::UTILS::matrix_rotate(vec, rot);
		ptr->set_coords(vec);
	}
}



void rotate_dihedral_forwards_ca(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives){

	//const int resCount = protein->get_residue_count();
	const int originResNum = angNum + 1;

	const vector3d origin = protein->get_bb_atom_coords(CA, originResNum);
	const PRODART::POSE::const_chain_shared_ptr orig_chain = protein->get_residue(originResNum)->get_chain();
	const int final_res_num = orig_chain->get_last_internal_residue_index();

	const vector3d origin_p1 = protein->get_bb_atom_coords(CA, originResNum + 1);


	const vector3d rot_axis = (origin_p1 - origin);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);

	for (int i = originResNum+1; i <= final_res_num; i++){
		atom_shared_ptr atm = protein->get_bb_atom(CA,i);
		if ((!run_through_inactives) && atm->isActive() == false) return;
		vector3d tempVec = atm->get_coords();//protein->get_bb_atom_coords(CA, i);
		//const PRODART::POSE::const_chain_shared_ptr curr_chain = protein->get_residue(i)->get_chain();
		//if (curr_chain != orig_chain) return;
		tempVec = tempVec - origin;
		quaternion_rotate(tempVec, params);
		tempVec = tempVec + origin;
		protein->set_bb_atom_coords(tempVec, CA, i);

	}

}

void rotate_dihedral_backwards_ca(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives){

	//const int resCount = protein->get_residue_count();
	const int originResNum = angNum + 2;

	const vector3d origin = protein->get_bb_atom_coords(CA, originResNum);
	const PRODART::POSE::const_chain_shared_ptr orig_chain = protein->get_residue(originResNum)->get_chain();
	const int first_res_num = orig_chain->get_first_internal_residue_index();

	const vector3d origin_m1 = protein->get_bb_atom_coords(CA, originResNum-1);

	const vector3d rot_axis = (origin_m1 - origin);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);

	for (int i = originResNum -1; i >= first_res_num; i--){
		atom_shared_ptr atm = protein->get_bb_atom(CA,i);
		if ((!run_through_inactives) && atm->isActive() == false) return;
		vector3d tempVec = atm->get_coords();//protein->get_bb_atom_coords(CA, i);
		//const PRODART::POSE::const_chain_shared_ptr curr_chain = protein->get_residue(i)->get_chain();
		//if (curr_chain != orig_chain) return;
		tempVec = tempVec - origin;
		quaternion_rotate(tempVec, params);
		tempVec = tempVec + origin;
		protein->set_bb_atom_coords(tempVec, CA, i);
	}


}


double rotate_dihedral_forwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives){

	//const int resCount = protein->get_residue_count();
	const int originResNum = angNum + 1;

	double max_move = 0;

	const vector3d origin = protein->get_bb_atom_coords(CA, originResNum);
	const PRODART::POSE::const_chain_shared_ptr orig_chain = protein->get_residue(originResNum)->get_chain();
	const int final_res_num = orig_chain->get_last_internal_residue_index();

	const vector3d origin_p1 = protein->get_bb_atom_coords(CA, originResNum + 1);


	const vector3d rot_axis = (origin_p1 - origin);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);

	for (int i = originResNum+1; i <= final_res_num; i++){
		atom_shared_ptr atm = protein->get_bb_atom(CA,i);
		if ((!run_through_inactives) && atm->isActive() == false) return std::sqrt(max_move);
		const vector3d startVec = atm->get_coords();//protein->get_bb_atom_coords(CA, i);
		//const PRODART::POSE::const_chain_shared_ptr curr_chain = protein->get_residue(i)->get_chain();
		//if (curr_chain != orig_chain) return;
		vector3d endVec = startVec - origin;
		quaternion_rotate(endVec, params);
		endVec = endVec + origin;
		protein->set_bb_atom_coords(endVec, CA, i);

		const double move_dist = (startVec - endVec).mod_sq();
		if (move_dist > max_move){
			max_move = move_dist;
		}

	}
	return std::sqrt(max_move);
}

double rotate_dihedral_backwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives){

	//const int resCount = protein->get_residue_count();
	const int originResNum = angNum + 2;

	double max_move = 0;

	const vector3d origin = protein->get_bb_atom_coords(CA, originResNum);
	const PRODART::POSE::const_chain_shared_ptr orig_chain = protein->get_residue(originResNum)->get_chain();
	const int first_res_num = orig_chain->get_first_internal_residue_index();

	const vector3d origin_m1 = protein->get_bb_atom_coords(CA, originResNum-1);

	const vector3d rot_axis = (origin_m1 - origin);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);

	for (int i = originResNum -1; i >= first_res_num; i--){
		atom_shared_ptr atm = protein->get_bb_atom(CA,i);
		if ((!run_through_inactives) && atm->isActive() == false) return std::sqrt(max_move);
		const vector3d startVec = atm->get_coords();//protein->get_bb_atom_coords(CA, i);
		//const PRODART::POSE::const_chain_shared_ptr curr_chain = protein->get_residue(i)->get_chain();
		//if (curr_chain != orig_chain) return;
		vector3d endVec = startVec - origin;
		quaternion_rotate(endVec, params);
		endVec = endVec + origin;
		protein->set_bb_atom_coords(endVec, CA, i);

		const double move_dist = (startVec - endVec).mod_sq();
		if (move_dist > max_move){
			max_move = move_dist;
		}

	}

	return std::sqrt(max_move);
}

double crankshaft_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int resNum1,
		const int resNum2,
		const double angle_rot,
		const bool run_through_inactives){


	const vector3d origin = protein->get_bb_atom_coords(CA, resNum1);//currentPdb->getVecByResidueIndex(resNum_i, CA, chain);

	const vector3d origin_p1 = protein->get_bb_atom_coords(CA, resNum2);//currentPdb->getVecByResidueIndex(resNum_j, CA, chain);


	const vector3d rot_axis = (origin_p1 - origin);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	double max_move = 0;

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);
	for (int i = resNum1 + 1; i < resNum2; i++){
		atom_shared_ptr atm = protein->get_bb_atom(CA,i);
		if ((!run_through_inactives) && atm->isActive() == false) return std::sqrt(max_move);
		const vector3d startVec = atm->get_coords();//currentPdb->getVecByResidueIndex(i, CA, chain);

		vector3d endVec = startVec - origin;
		quaternion_rotate(endVec, params);
		endVec = endVec + origin;
		atm->set_coords(endVec);
		//currentPdb->setVecByResidueIndex(i, CA, chain, endVec);
		const double move_dist = (startVec - endVec).mod_sq();
		if (move_dist > max_move){
			max_move = move_dist;
		}


	}




	return std::sqrt(max_move);
}


double rotate_angle_forwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives ){

	double max_move = 0;


	const int originResNum = angNum + 1;

	const vector3d origin =  protein->get_bb_atom_coords(CA, originResNum); // currentPdb->getVecByResidueIndex(originResNum, CA, chain);

	const vector3d origin_p1 = protein->get_bb_atom_coords(CA, originResNum+1); //currentPdb->getVecByResidueIndex(originResNum+1, CA, chain);

	const vector3d origin_m1 = protein->get_bb_atom_coords(CA, originResNum-1); //currentPdb->getVecByResidueIndex(originResNum-1, CA, chain);

	const vector3d rot_axis = (origin_p1 - origin) * (origin - origin_m1);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	const PRODART::POSE::const_chain_shared_ptr orig_chain = protein->get_residue(originResNum)->get_chain();
	const int final_res_num = orig_chain->get_last_internal_residue_index();

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);
	for (int i = originResNum; i <= final_res_num; i++){
		atom_shared_ptr atm = protein->get_bb_atom(CA,i);
		if ((!run_through_inactives) && atm->isActive() == false) return std::sqrt(max_move);
		const vector3d startVec = atm->get_coords();
		vector3d endVec = startVec - origin;
		quaternion_rotate(endVec, params);
		endVec = endVec + origin;
		atm->set_coords(endVec);
		const double move_dist = (startVec - endVec).mod_sq();
		if (move_dist > max_move){
			max_move = move_dist;
		}
	}


	return std::sqrt(max_move);//true;
}

double rotate_angle_backwards_ca_track_move(const PRODART::POSE::pose_shared_ptr protein,
		const int angNum,
		const double angle_rot,
		const bool run_through_inactives ){
	const int originResNum = angNum + 1;
	double max_move = 0;

	const vector3d origin = protein->get_bb_atom_coords(CA, originResNum); // currentPdb->getVecByResidueIndex(originResNum, CA, chain);

	const vector3d origin_p1 = protein->get_bb_atom_coords(CA, originResNum+1); //currentPdb->getVecByResidueIndex(originResNum+1, CA, chain);

	const vector3d origin_m1 = protein->get_bb_atom_coords(CA, originResNum-1); //currentPdb->getVecByResidueIndex(originResNum-1, CA, chain);


	const vector3d rot_axis = (origin_p1 - origin) * (origin - origin_m1);
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();

	const PRODART::POSE::const_chain_shared_ptr orig_chain = protein->get_residue(originResNum)->get_chain();
	const int first_res_num = orig_chain->get_first_internal_residue_index();


	//cout << rot_axis_norm << endl;

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);
	for (int i = originResNum -1; i >= first_res_num; i--){
		atom_shared_ptr atm = protein->get_bb_atom(CA,i);
		if ((!run_through_inactives) && atm->isActive() == false) return std::sqrt(max_move);
		const vector3d startVec = atm->get_coords();
		vector3d endVec = startVec - origin;
		quaternion_rotate(endVec, params);
		endVec = endVec + origin;

		atm->set_coords(endVec);
		const double move_dist = (startVec - endVec).mod_sq();
		if (move_dist > max_move){
			max_move = move_dist;
		}

	}

	return std::sqrt(max_move);
}

bool check_chi_def(const atom_shared_ptr_vector& vec){

	if (vec.size() != 4){
		return false;
	}

	atom_shared_ptr_vector::const_iterator it;
	for (it = vec.begin(); it != vec.end(); it++){
		if (!(*it)){
			return false;
		}
		if (!(*it)->isActiveAndSet()){
			return false;
		}
	}

	return true;
}


bool set_chi_angle_forwards(const PRODART::POSE::atom_shared_ptr_vector& chi_def,
		const PRODART::POSE::atom_shared_ptr_vector& chi_fwd,
		const double chi,
		const bool run_through_inactives){
	//boost::tuple<const atom_shared_ptr_vector&, const atom_shared_ptr_vector&, const atom_shared_ptr_vector&> details = sc->get_chi_details(chi_index);

	if (check_chi_def(chi_def) != true){
		return false;
	}

	const double old_chi = UTILS::dihedral(chi_def[0]->get_coords(), chi_def[1]->get_coords(), chi_def[2]->get_coords(), chi_def[3]->get_coords());//sc->get_chi(chi_index);
	const double angle_rot = chi - old_chi;

	const vector3d rot_axis = chi_def[2]->get_coords() - chi_def[1]->get_coords();
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();
	const vector3d origin = chi_def[2]->get_coords();

	/*
	for (int i = 0; i< 4; i++){
		cout << details.get<0>()[i]->get_type().get_label() << " ";
	}
	cout << endl;
	*/

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);
	atom_shared_ptr_vector::const_iterator it;
	for (it = chi_fwd.begin(); it != chi_fwd.end(); it++){
		atom_shared_ptr atm = *it;
		if (atm && (atm->isActiveAndSet() || run_through_inactives)){
			const vector3d startVec = atm->get_coords();
			vector3d endVec = startVec - origin;
			quaternion_rotate(endVec, params);
			endVec = endVec + origin;
			atm->set_coords(endVec);
			//cout << atm->get_type().get_label() << " ";
		}


	}
	//cout << endl;

#ifndef NDEBUG

	const double new_chi = UTILS::dihedral(chi_def[0]->get_coords(), chi_def[1]->get_coords(), chi_def[2]->get_coords(), chi_def[3]->get_coords());//sc->get_chi(chi_index);
	if (fabs(UTILS::get_dihedral_distance(new_chi, chi)) > 10e-2 ){
		std::cerr << "set_chi_angle_forwards: ERROR:\t" << UTILS::radians_to_degrees(old_chi) << "\t"<< UTILS::radians_to_degrees(chi) << "\t" << UTILS::radians_to_degrees(new_chi) << endl;
	}
	assert( fabs(UTILS::get_dihedral_distance(new_chi, chi)) < 10e-2  );

#endif


	return true;
}

bool set_chi_angle_backwards(const PRODART::POSE::atom_shared_ptr_vector& chi_def,
		const PRODART::POSE::atom_shared_ptr_vector& chi_bwd,
		const double chi,
		const bool run_through_inactives){
	//boost::tuple<const atom_shared_ptr_vector&, const atom_shared_ptr_vector&, const atom_shared_ptr_vector&> details = sc->get_chi_details(chi_index);

	if (check_chi_def(chi_def) != true){
		return false;
	}

	const double old_chi = UTILS::dihedral(chi_def[0]->get_coords(), chi_def[1]->get_coords(), chi_def[2]->get_coords(), chi_def[3]->get_coords());//sc->get_chi(chi_index);
	const double angle_rot = chi - old_chi;

	const vector3d rot_axis = chi_def[1]->get_coords() - chi_def[2]->get_coords();
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();
	const vector3d origin = chi_def[1]->get_coords();

	/*
	for (int i = 0; i< 4; i++){
		cout << details.get<0>()[i]->get_type().get_label() << " ";
	}
	cout << endl;
	*/

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);
	atom_shared_ptr_vector::const_iterator it;
	for (it = chi_bwd.begin(); it != chi_bwd.end(); it++){
		atom_shared_ptr atm = *it;
		if (atm && (atm->isActiveAndSet() || run_through_inactives)){
			const vector3d startVec = atm->get_coords();
			vector3d endVec = startVec - origin;
			quaternion_rotate(endVec, params);
			endVec = endVec + origin;
			atm->set_coords(endVec);
			//cout << atm->get_type().get_label() << " ";
		}
		/*
		if (!atm->isActiveAndSet()){
			cout << "atom " << atm->get_type().get_label() << " is inactive" << endl;
		}
		*/

	}
	//cout << endl;

#ifndef NDEBUG

	const double new_chi = UTILS::dihedral(chi_def[0]->get_coords(), chi_def[1]->get_coords(), chi_def[2]->get_coords(), chi_def[3]->get_coords());//sc->get_chi(chi_index);
	if (fabs(UTILS::get_dihedral_distance(new_chi, chi)) > 10e-2 ){
		std::cerr << "set_chi_angle_backwards: ERROR:\t" << UTILS::radians_to_degrees(old_chi) << "\t"<< UTILS::radians_to_degrees(chi) << "\t" << UTILS::radians_to_degrees(new_chi) << endl;
	}
	assert( fabs(UTILS::get_dihedral_distance(new_chi, chi)) < 10e-2  );

#endif




	return true;
}

bool set_chi_angle_forwards(PRODART::POSE::sidechain_shared_ptr sc,
		const unsigned int chi_index,
		const double chi,
		const bool run_through_inactives){

	boost::tuple<const atom_shared_ptr_vector&, const atom_shared_ptr_vector&, const atom_shared_ptr_vector&> details = sc->get_chi_details(chi_index);


	/*
	if (check_chi_def(details.get<0>()) != true){
		return false;
	}

	const double old_chi = sc->get_chi(chi_index);
	const double angle_rot = chi - old_chi;

	const vector3d rot_axis = details.get<0>()[2]->get_coords() - details.get<0>()[1]->get_coords();
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();
	const vector3d origin = details.get<0>()[2]->get_coords();

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);
	atom_shared_ptr_vector::const_iterator it;
	for (it = details.get<1>().begin(); it != details.get<1>().end(); it++){
		atom_shared_ptr atm = *it;
		if (atm && (atm->isActiveAndSet() || run_through_inactives)){
			const vector3d startVec = atm->get_coords();
			vector3d endVec = startVec - origin;
			quaternion_rotate(endVec, params);
			endVec = endVec + origin;
			atm->set_coords(endVec);
			//cout << atm->get_type().get_label() << " ";
		}


	}

	//cout << endl;

#ifndef NDEBUG

	const double new_chi = sc->get_chi(chi_index);
	if (fabs(UTILS::get_dihedral_distance(new_chi, chi)) > 10e-2 ){
		std::cerr << "set_chi_angle_forwards: ERROR:\t" << UTILS::radians_to_degrees(old_chi) << "\t"<< UTILS::radians_to_degrees(chi) << "\t" << UTILS::radians_to_degrees(new_chi) << endl;
	}
	assert( fabs(UTILS::get_dihedral_distance(new_chi, chi)) < 10e-2  );

#endif
	 */

	return set_chi_angle_forwards(details.get<0>(), details.get<1>(), chi, run_through_inactives);

	//return true;
}

bool set_chi_angle_backwards(PRODART::POSE::sidechain_shared_ptr sc,
		const unsigned int chi_index,
		const double chi,
		const bool run_through_inactives){
	boost::tuple<const atom_shared_ptr_vector&, const atom_shared_ptr_vector&, const atom_shared_ptr_vector&> details = sc->get_chi_details(chi_index);


	/*
	if (check_chi_def(details.get<0>()) != true){
		return false;
	}

	const double old_chi = sc->get_chi(chi_index);
	const double angle_rot = chi - old_chi;

	const vector3d rot_axis = details.get<0>()[1]->get_coords() - details.get<0>()[2]->get_coords();
	const vector3d rot_axis_norm = rot_axis / rot_axis.mod();
	const vector3d origin = details.get<0>()[1]->get_coords();

	const quaternion_rotate_params params(angle_rot, rot_axis_norm);
	atom_shared_ptr_vector::const_iterator it;
	for (it = details.get<2>().begin(); it != details.get<2>().end(); it++){
		atom_shared_ptr atm = *it;
		if (atm && (atm->isActiveAndSet() || run_through_inactives)){
			const vector3d startVec = atm->get_coords();
			vector3d endVec = startVec - origin;
			quaternion_rotate(endVec, params);
			endVec = endVec + origin;
			atm->set_coords(endVec);
			//cout << atm->get_type().get_label() << " ";
		}


	}
	//cout << endl;

#ifndef NDEBUG

	const double new_chi = sc->get_chi(chi_index);
	if (fabs(UTILS::get_dihedral_distance(new_chi, chi)) > 10e-2 ){
		std::cerr << "set_chi_angle_backwards: ERROR:\t" << UTILS::radians_to_degrees(old_chi) << "\t"<< UTILS::radians_to_degrees(chi) << "\t" << UTILS::radians_to_degrees(new_chi) << endl;
	}
	assert( fabs(UTILS::get_dihedral_distance(new_chi, chi)) < 10e-2  );

#endif
	 */

	return set_chi_angle_backwards(details.get<0>(), details.get<2>(), chi, run_through_inactives);


	return true;
}


}
}
}
