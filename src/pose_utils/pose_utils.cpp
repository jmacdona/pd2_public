//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose_utils.cpp
 *
 *  Created on: Feb 19, 2010
 *      Author: jmacdona
 */

#include "pose_utils.h"

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

namespace PRODART {
namespace POSE_UTILS {

typedef std::vector<string> string_vector;

using namespace KINE;

bool get_line_CA_segment(const PRODART::POSE::const_pose_shared_ptr protein,
		const int start_residue_index,
		const int end_residue_index,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec){

	const int seg_len = end_residue_index - start_residue_index + 1;

	if (seg_len < 1){
		return false;
	}
	else if (seg_len == 1){
		startVec = protein->get_bb_atom_coords(CA, start_residue_index);
		endVec = protein->get_bb_atom_coords(CA, start_residue_index);
		return false;
	}
	else if (seg_len == 2){
		startVec = protein->get_bb_atom_coords(CA, start_residue_index);
		endVec = protein->get_bb_atom_coords(CA, end_residue_index);
		return true;
	}

	PRODART::UTILS::vector3d_vector coords;
	coords.reserve(seg_len);

	for (int i = start_residue_index; i <= end_residue_index; i++){

		const vector3d ca_vec = protein->get_bb_atom_coords(CA, i);

		coords.push_back(ca_vec);


	}

	PRODART::UTILS::line_fit3d(coords,
			startVec,
			endVec);


	return true;

}




double get_ca_rmsd(const PRODART::POSE::const_pose_shared_ptr protein1,
		const PRODART::POSE::const_pose_shared_ptr protein2){

	const int resCount1 = protein1->get_residue_count();
	const int resCount2 = protein2->get_residue_count();

	if (resCount1 != resCount2){
		std::cerr << "ERROR: get_ca_rmsd can not calculate RMSD between different size proteins" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	double **pdbA_coords = MatInit(3, resCount1), **pdbB_coords = MatInit(3, resCount1);
    //double rotmat[9];
	for(int i = 0; i < resCount1; ++ i){
		//const int i_pos =  (i * PRODART::PROT::num_mainchain_atoms_per_residue) + CA_pos;
		/** pdb A */
		//if (protein1->get_bb_atom(CA, i)->isSet() && protein2->get_bb_atom(CA, i)->isSet()){
		pdbA_coords[0][i] =  protein1->get_bb_atom_coords(CA, i)[0];
		pdbA_coords[1][i] =  protein1->get_bb_atom_coords(CA, i)[1];
		pdbA_coords[2][i] =  protein1->get_bb_atom_coords(CA, i)[2];
		/** pdb B */
		pdbB_coords[0][i] =  protein2->get_bb_atom_coords(CA, i)[0];
		pdbB_coords[1][i] =  protein2->get_bb_atom_coords(CA, i)[1];
		pdbB_coords[2][i] =  protein2->get_bb_atom_coords(CA, i)[2];
		//}
	}


	const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, resCount1);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	return rmsd;


}

double get_ca_rmsd_superpose(const PRODART::POSE::const_pose_shared_ptr protein1,
		const PRODART::POSE::pose_shared_ptr protein2){
	const int resCount1 = protein1->get_residue_count();
	const int resCount2 = protein2->get_residue_count();

	if (resCount1 != resCount2){
		std::cerr << "ERROR: get_ca_rmsd can not calculate RMSD between different size proteins" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	double **pdbA_coords = MatInit(3, resCount1), **pdbB_coords = MatInit(3, resCount1);
    double rotmat[9];
	for(int i = 0; i < resCount1; ++ i){
		//const int i_pos =  (i * PRODART::PROT::num_mainchain_atoms_per_residue) + CA_pos;
		/** pdb A */
		//if (protein1->get_bb_atom(CA, i)->isSet() && protein2->get_bb_atom(CA, i)->isSet()){
		pdbA_coords[0][i] =  protein1->get_bb_atom_coords(CA, i)[0];
		pdbA_coords[1][i] =  protein1->get_bb_atom_coords(CA, i)[1];
		pdbA_coords[2][i] =  protein1->get_bb_atom_coords(CA, i)[2];
		/** pdb B */
		pdbB_coords[0][i] =  protein2->get_bb_atom_coords(CA, i)[0];
		pdbB_coords[1][i] =  protein2->get_bb_atom_coords(CA, i)[1];
		pdbB_coords[2][i] =  protein2->get_bb_atom_coords(CA, i)[2];
		//}
	}

	vector3d prot1_com ;//= get_ca_CoM(protein1);
	vector3d prot2_com ;//= get_ca_CoM(protein2);

	const double rmsd = CalcRMSDRotationalMatrix(pdbA_coords, pdbB_coords, resCount1,
			prot1_com, prot2_com,
			rotmat, NULL);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);



	PRODART::UTILS::rot_matrix rot(rotmat);

	translate_fa(protein2, -prot2_com);
	apply_rotation_matrix_fa(protein2, rot);
	translate_fa(protein2, prot1_com);

	return rmsd;
}

double get_ca_rmsd(const PRODART::POSE::const_pose_shared_ptr protein1,
		const PRODART::POSE::const_pose_shared_ptr protein2,
		const std::map<int,int> residue_mapping){

	const int aln_size = residue_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_ca_rmsd can not calculate RMSD without an alignment" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);

	std::map<int,int>::const_iterator iter;
	int i = 0;
	for (iter = residue_mapping.begin(); iter != residue_mapping.end(); iter++){
		// pdb A /
		pdbA_coords[0][i] =  protein1->get_bb_atom_coords(CA, iter->first)[0];
		pdbA_coords[1][i] =  protein1->get_bb_atom_coords(CA, iter->first)[1];
		pdbA_coords[2][i] =  protein1->get_bb_atom_coords(CA, iter->first)[2];
		// pdb B /
		pdbB_coords[0][i] =  protein2->get_bb_atom_coords(CA, iter->second)[0];
		pdbB_coords[1][i] =  protein2->get_bb_atom_coords(CA, iter->second)[1];
		pdbB_coords[2][i] =  protein2->get_bb_atom_coords(CA, iter->second)[2];
		i++;
	}

	const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);



	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);



	return rmsd;

}

double get_ca_rmsd_superpose(const PRODART::POSE::const_pose_shared_ptr protein1,
		const PRODART::POSE::pose_shared_ptr protein2,
		const std::map<int,int> residue_mapping){

	const int aln_size = residue_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_ca_rmsd can not calculate RMSD without an alignment" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);
    double rotmat[9];
	std::map<int,int>::const_iterator iter;
	int i = 0;
	for (iter = residue_mapping.begin(); iter != residue_mapping.end(); iter++){
		// pdb A /
		pdbA_coords[0][i] =  protein1->get_bb_atom_coords(CA, iter->first)[0];
		pdbA_coords[1][i] =  protein1->get_bb_atom_coords(CA, iter->first)[1];
		pdbA_coords[2][i] =  protein1->get_bb_atom_coords(CA, iter->first)[2];
		// pdb B /
		pdbB_coords[0][i] =  protein2->get_bb_atom_coords(CA, iter->second)[0];
		pdbB_coords[1][i] =  protein2->get_bb_atom_coords(CA, iter->second)[1];
		pdbB_coords[2][i] =  protein2->get_bb_atom_coords(CA, iter->second)[2];
		i++;
	}

	vector3d prot1_com ;//= get_ca_CoM(protein1);
	vector3d prot2_com ;//= get_ca_CoM(protein2);

	const double rmsd = CalcRMSDRotationalMatrix(pdbA_coords, pdbB_coords, aln_size,
			prot1_com, prot2_com,
			rotmat, NULL);
	//const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);
	PRODART::UTILS::rot_matrix rot(rotmat);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	translate_fa(protein2, -prot2_com);
	apply_rotation_matrix_fa(protein2, rot);
	translate_fa(protein2, prot1_com);

	return rmsd;

}

double get_rmsd(const std::map<PRODART::POSE::const_atom_shared_ptr,PRODART::POSE::const_atom_shared_ptr> atom_mapping){

	const int aln_size = atom_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_rmsd can not calculate RMSD without an alignment" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);

	std::map<const_atom_shared_ptr,const_atom_shared_ptr>::const_iterator iter;
	int i = 0;
	for (iter = atom_mapping.begin(); iter != atom_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->first)->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->second)->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		i++;
	}

	const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	return rmsd;

}

double get_msd(const std::map<PRODART::POSE::const_atom_shared_ptr,PRODART::POSE::const_atom_shared_ptr> atom_mapping){

	const int aln_size = atom_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_msd can not calculate RMSD without an alignment" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);

	std::map<const_atom_shared_ptr,const_atom_shared_ptr>::const_iterator iter;
	int i = 0;
	for (iter = atom_mapping.begin(); iter != atom_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->first)->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->second)->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		i++;
	}

	const double msd = QCP_msd(pdbA_coords, pdbB_coords, aln_size);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	return msd;

}

double get_msd(const std::vector< boost::tuple <PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr, double > > atom_atom_wt_mapping,
		const PRODART::UTILS::vector3d& CoM1, const PRODART::UTILS::vector3d& CoM2){
	const int aln_size = atom_atom_wt_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_msd can not calculate RMSD without an alignment" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}
	/*
	else if (aln_size != (int)weights.size()){
		std::cerr << "ERROR: get_msd: weights not same size as atom_mapping" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	*/

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);

	std::vector< boost::tuple <PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr, double > >::const_iterator iter;
	double cweights[aln_size];
	int i = 0;
	for (iter = atom_atom_wt_mapping.begin(); iter != atom_atom_wt_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->get<0>())->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->get<1>())->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		cweights[i] = (iter->get<2>());
		//cout << "DEBUGGER-get_msd\t" << atomA << atomB << endl;
		i++;
	}


	const double msd = QCP_msd_force_CoM(pdbA_coords, pdbB_coords, aln_size, CoM1, CoM2, cweights);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	return msd;
}


double get_rmsd(const std::map<PRODART::POSE::atom_shared_ptr,PRODART::POSE::atom_shared_ptr> atom_mapping){

	const int aln_size = atom_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_rmsd can not calculate RMSD without an alignment" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);

	std::map<atom_shared_ptr,atom_shared_ptr>::const_iterator iter;
	int i = 0;
	for (iter = atom_mapping.begin(); iter != atom_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->first)->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->second)->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		i++;
	}

	const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	return rmsd;

}

double get_msd(const std::map<PRODART::POSE::atom_shared_ptr,PRODART::POSE::atom_shared_ptr> atom_mapping){

	const int aln_size = atom_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_msd can not calculate RMSD without an alignment" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);

	std::map<atom_shared_ptr,atom_shared_ptr>::const_iterator iter;
	int i = 0;
	for (iter = atom_mapping.begin(); iter != atom_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->first)->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->second)->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		i++;
	}

	const double msd = QCP_msd(pdbA_coords, pdbB_coords, aln_size);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	return msd;

}

double get_rmsd_superpose(const PRODART::POSE::pose_shared_ptr protein_to_move,
		const std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping){
	const int aln_size = atom_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_rmsd_superpose can not calculate RMSD without an atom mapping" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);
    double rotmat[9];
	std::map<const_atom_shared_ptr,const_atom_shared_ptr>::const_iterator iter;
	int i = 0;
	for (iter = atom_mapping.begin(); iter != atom_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->first)->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->second)->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		i++;
	}

	vector3d prot1_com ;//= get_ca_CoM(protein1);
	vector3d prot2_com ;//= get_ca_CoM(protein2);

	const double rmsd = CalcRMSDRotationalMatrix(pdbA_coords, pdbB_coords, aln_size,
			prot1_com, prot2_com,
			rotmat, NULL);
	//const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);
	PRODART::UTILS::rot_matrix rot(rotmat);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	translate_fa(protein_to_move, -prot2_com);
	apply_rotation_matrix_fa(protein_to_move, rot);
	translate_fa(protein_to_move, prot1_com);

	return rmsd;
}

double get_rmsd_superpose(const PRODART::POSE::pose_shared_ptr protein_to_move,
		const std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping){
	const int aln_size = atom_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_rmsd_superpose can not calculate RMSD without an atom mapping" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);
    double rotmat[9];
	std::map<atom_shared_ptr,atom_shared_ptr>::const_iterator iter;
	int i = 0;
	for (iter = atom_mapping.begin(); iter != atom_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->first)->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->second)->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		i++;
	}

	vector3d prot1_com ;//= get_ca_CoM(protein1);
	vector3d prot2_com ;//= get_ca_CoM(protein2);

	const double rmsd = CalcRMSDRotationalMatrix(pdbA_coords, pdbB_coords, aln_size,
			prot1_com, prot2_com,
			rotmat, NULL);
	//const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);
	PRODART::UTILS::rot_matrix rot(rotmat);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	translate_fa(protein_to_move, -prot2_com);
	apply_rotation_matrix_fa(protein_to_move, rot);
	translate_fa(protein_to_move, prot1_com);

	return rmsd;
}

double get_rmsd_superpose(const PRODART::POSE::pose_shared_ptr protein_to_move,
		const std::vector<boost::tuple<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr, double> > atom_atom_wt_mapping,
		const PRODART::UTILS::vector3d& CoM_ref, const PRODART::UTILS::vector3d& CoM_to_move){
	const int aln_size = atom_atom_wt_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_rmsd_superpose can not calculate RMSD without an atom mapping" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}
	/*
	else if (aln_size != (int)weights.size()){
		std::cerr << "ERROR: get_rmsd_superpose: weights not same size as atom_mapping" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	*/

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);
    double rotmat[9];
    std::vector<boost::tuple<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr, double> >::const_iterator iter;
	double cweights[atom_atom_wt_mapping.size()];
	int i = 0;
	for (iter = atom_atom_wt_mapping.begin(); iter != atom_atom_wt_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->get<0>())->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->get<1>())->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		cweights[i] = (iter->get<2>());
		i++;
	}

	//vector3d prot1_com ;//= get_ca_CoM(protein1);
	//vector3d prot2_com ;//= get_ca_CoM(protein2);

/*
	for (unsigned int i = 0; i < weights.size(); i++){
		cweights[i] = weights[i];
	}
	*/

	const double rmsd = CalcRMSDRotationalMatrix_force_CoM(pdbA_coords, pdbB_coords, aln_size,
			CoM_ref, CoM_to_move,
			rotmat, cweights);
	//const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);
	PRODART::UTILS::rot_matrix rot(rotmat);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	translate_fa(protein_to_move, -CoM_to_move);
	apply_rotation_matrix_fa(protein_to_move, rot);
	translate_fa(protein_to_move, CoM_ref);

	return rmsd;
}

//! calculate RMSD given an arbitrary atom_mapping - can be used with any atoms - then translates and rotates the supplied atoms (should but doesn't have to correspond to 'second' in the atom map) - NOT rigorously tested yet
double get_rmsd_superpose(const PRODART::POSE::atom_shared_ptr_vector atoms_to_move,
		const std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> atom_mapping){
	const int aln_size = atom_mapping.size();

	if (aln_size == 0){
		std::cerr << "ERROR: get_rmsd_superpose can not calculate RMSD without an atom mapping" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if(aln_size == 1){
		return 0;
	}

	double **pdbA_coords = MatInit(3, aln_size), **pdbB_coords = MatInit(3, aln_size);
    double rotmat[9];
	std::map<atom_shared_ptr,atom_shared_ptr>::const_iterator iter;
	int i = 0;
	for (iter = atom_mapping.begin(); iter != atom_mapping.end(); iter++){
		// pdb A /
		const vector3d atomA = (iter->first)->get_coords();
		pdbA_coords[0][i] =  atomA[0];
		pdbA_coords[1][i] =  atomA[1];
		pdbA_coords[2][i] =  atomA[2];
		// pdb B /
		const vector3d atomB = (iter->second)->get_coords();
		pdbB_coords[0][i] =  atomB[0];
		pdbB_coords[1][i] =  atomB[1];
		pdbB_coords[2][i] =  atomB[2];
		i++;
	}

	vector3d prot1_com ;//= get_ca_CoM(protein1);
	vector3d prot2_com ;//= get_ca_CoM(protein2);

	const double rmsd = CalcRMSDRotationalMatrix(pdbA_coords, pdbB_coords, aln_size,
			prot1_com, prot2_com,
			rotmat, NULL);
	//const double rmsd = QCP_rmsd(pdbA_coords, pdbB_coords, aln_size);
	PRODART::UTILS::rot_matrix rot(rotmat);

	MatDestroy(&pdbA_coords);
	MatDestroy(&pdbB_coords);

	/*
	translate_fa(protein_to_move, -prot2_com);
	apply_rotation_matrix_fa(protein_to_move, rot);
	translate_fa(protein_to_move, prot1_com);
	*/
	// apply translation and rotation to  atoms_to_move instead
	translate_fa(atoms_to_move, -prot2_com);
	apply_rotation_matrix_fa(atoms_to_move, rot);
	translate_fa(atoms_to_move, prot1_com);

	return rmsd;
}










double get_fixed_position_rmsd_auto_loop(const PRODART::POSE::const_pose_shared_ptr pose1,
		const PRODART::POSE::pose_shared_ptr pose2,
		const double tol){

	const int rescount1 = pose1->get_residue_count();
	const int rescount2 = pose2->get_residue_count();

	if (rescount1 != rescount2){
		cerr << "get_fixed_position_rmsd_auto_loop: ERROR: poses are not identical in length" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}


	std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> atom_mapping;

	for (int resnum = 0; resnum < rescount1; resnum++){
		const int start_bb_at = pose1->get_residue(resnum)->get_first_bb_atom_index();
		for (int bb_atom_num = start_bb_at; bb_atom_num < start_bb_at + pose::get_num_protein_backbone_atom_types(); bb_atom_num++ ){
			const_atom_shared_ptr ref_atom = pose1->get_bb_atom(bb_atom_num);
			const_atom_shared_ptr tm_atom = pose2->get_atom(ref_atom->get_type(),resnum);
			if (ref_atom && tm_atom){
				if (ref_atom->isSet() && tm_atom->isSet()){
					if ((ref_atom->get_coords() - tm_atom->get_coords()).mod() >= tol){
						atom_mapping[ref_atom] = tm_atom;
						cout << "mapped:\t"
								<< ref_atom->get_residue()->get_type().get_label() << "\t"
								<< ref_atom->get_residue()->get_pdb_residue_index() << "\t"
								<< ref_atom->get_chain()->getChainID() << "\t"
								<< ref_atom->get_type().get_pdb_formatted_label() << "\t"
								<< "->" << "\t"
								<< tm_atom->get_residue()->get_type().get_label() << "\t"
								<< tm_atom->get_residue()->get_pdb_residue_index() << "\t"
								<< tm_atom->get_chain()->getChainID() << "\t"
								<< tm_atom->get_type().get_pdb_formatted_label()
								<< endl;
					}
				}

			}
		}
		//  add in sidechain atom mapping here
		const int num_sc_atoms = pose1->get_residue(resnum)->get_sidechain()->get_atom_count();
		for (int i = 0; i < num_sc_atoms; i++){
			const_atom_shared_ptr ref_atom = pose1->get_residue(resnum)->get_sidechain()->get_atom(i);
			const_atom_shared_ptr tm_atom = pose2->get_atom(ref_atom->get_type(),resnum);
			if (ref_atom && tm_atom){
				if (ref_atom->isSet() && tm_atom->isSet()){
					if ((ref_atom->get_coords() - tm_atom->get_coords()).mod() >= tol){
						atom_mapping[ref_atom] = tm_atom;
						cout << "mapped:\t"
								<< ref_atom->get_residue()->get_type().get_label() << "\t"
								<< ref_atom->get_residue()->get_pdb_residue_index() << "\t"
								<< ref_atom->get_chain()->getChainID() << "\t"
								<< ref_atom->get_type().get_pdb_formatted_label() << "\t"
								<< "->" << "\t"
								<< tm_atom->get_residue()->get_type().get_label() << "\t"
								<< tm_atom->get_residue()->get_pdb_residue_index() << "\t"
								<< tm_atom->get_chain()->getChainID() << "\t"
								<< tm_atom->get_type().get_pdb_formatted_label()
								<< endl;
					}
				}

			}

		}
	}

	cout << "total atoms mapped:\t"  << atom_mapping.size() << endl;

	return get_fixed_position_rmsd(atom_mapping);
}

double get_fixed_position_rmsd(const std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr> &atom_mapping){
	if (atom_mapping.size() != 0){
		std::map<PRODART::POSE::const_atom_shared_ptr, PRODART::POSE::const_atom_shared_ptr>::const_iterator it;
		double sum_sq = 0;
		for (it = atom_mapping.begin(); it!=atom_mapping.end(); it++){
			sum_sq += (it->first->get_coords() - it->second->get_coords()).mod_sq();
		}
		const double len =  (double)atom_mapping.size();
		return sqrt(sum_sq / len);
	}
	else {
		return std::numeric_limits<double>::quiet_NaN();
	}
}

double get_fixed_position_rmsd(PRODART::POSE::atom_shared_ptr_vector& vec1, PRODART::POSE::atom_shared_ptr_vector& vec2){
	if (vec1.size() != vec2.size()){
		cerr << "get_fixed_position_rmsd: ERROR: atoms vectors are not identical in length" << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	double sum_sq = 0;
	for (unsigned int i = 0; i < vec1.size(); i++){
		sum_sq += (vec1[i]->get_coords() - vec2[i]->get_coords()).mod_sq();
	}
	const double len =  (double)vec1.size();

	return sqrt(sum_sq / len);
}

//TODO deal with isSetActive non proteins
PRODART::UTILS::vector3d get_ca_CoM(const PRODART::POSE::const_pose_shared_ptr protein){
	vector3d vec_sum(0,0,0);

	const int resCount = protein->get_residue_count();
	double real_count = 0;

	for (int i = 0; i < resCount; i++){
		const_atom_shared_ptr atm = protein->get_bb_atom(CA, i);
		if (atm->isActiveAndSet()){
			const_chain_shared_ptr ch = atm->get_chain();
			if (ch->isPeptide()){
				const vector3d vec = protein->get_bb_atom_coords(CA, i);
				real_count += 1.0;
				vec_sum += vec;
			}
		}

	}
	const double numRes = real_count;// static_cast<double>(resCount);

	vector3d mean_vec = vec_sum / numRes;
	return mean_vec;
}


double get_ca_radgyr(const PRODART::POSE::const_pose_shared_ptr protein,
		const int startRes,
		const int endRes){
	vector3d vec_sum(0,0,0);
	vector3d vec_sum_sq(0,0,0);

	for (int i = startRes; i <= endRes; i++){
		const vector3d vec = protein->get_bb_atom_coords(CA, i);

		vec_sum += vec;
		vec_sum_sq += vec.sq();
	}
	const double numRes = static_cast<double>(endRes - startRes + 1);

	vector3d mean_vec = vec_sum / numRes;
	vector3d var_vec = (vec_sum_sq / numRes) - mean_vec.sq();

	//sum = sum / (numRes * numRes);
	return std::sqrt(var_vec[0] + var_vec[1] + var_vec[2]);
}

double get_ca_radgyr(const PRODART::POSE::const_pose_shared_ptr protein,
		const int chain_num){
	const_chain_shared_ptr currChain = protein->get_chain(chain_num);
	const int startRes = currChain->get_first_internal_residue_index();
	const int endRes = currChain->get_last_internal_residue_index();

	return get_ca_radgyr(protein, startRes, endRes);
}


//! quick and dirty addition of backbone hydrogen to a residue
UTILS::vector3d get_ideal_HN(const PRODART::POSE::pose_shared_ptr protein,
		const int residue_number,
		bool &is_defined){

	is_defined = true;
	//cout << "called quick_add_HN for res num:" << residue_number << "\t";

	const residue_shared_ptr currRes = protein->get_residue(residue_number);
	const const_residue_shared_ptr prevRes = currRes->get_prev_residue();



	if (currRes && prevRes && residue_number > 0){
		const_atom_shared_ptr c_m1_ptr = protein->get_bb_atom(C, residue_number - 1),
				ca_0_ptr = protein->get_bb_atom(CA, residue_number),
				n_0_ptr = protein->get_bb_atom(N, residue_number);
		if (c_m1_ptr->isSet() && ca_0_ptr->isSet() && n_0_ptr->isSet()){
			const vector3d new_h_coords = PRODART::UTILS::dihedralEnd(c_m1_ptr->get_coords(),
					ca_0_ptr->get_coords(),
					n_0_ptr->get_coords(),
					0.9992, 								// bond_3_4
					115.4200 * (PRODART::UTILS::PI / 180.0),				// angle_2_3_4
					180.0000 * (PRODART::UTILS::PI / 180.0));				// dihedral_1_2_3_4

			//cout << "!!!" << h_0_ptr->get_coords() << "!!!" << endl;
			return new_h_coords;
		}
		else {
			// not enough atoms - figure out something else here
			is_defined = false;

			//cout << "not enough atoms\n";


		}
	}
	else {
		// N-terminal residue - figure out something else here
		// undefined
		const_atom_shared_ptr c_0_ptr = protein->get_bb_atom(C, residue_number),
				ca_0_ptr = protein->get_bb_atom(CA, residue_number),
				n_0_ptr = protein->get_bb_atom(N, residue_number);

		if ( !c_0_ptr->isSet() || !ca_0_ptr->isSet()
				|| !n_0_ptr->isSet() ) {
			is_defined = false;
			return vector3d();
		}

		double angle = 0;

		//srand(10);
		angle = ((180.0/180.0) * PRODART::UTILS::PI);


		const vector3d new_h_coords = dihedralEnd( c_0_ptr->get_coords(),
				ca_0_ptr->get_coords(),
				n_0_ptr->get_coords(),
				0.9992,
				((115.0/180.0) * PRODART::UTILS::PI) ,
				angle );

		return new_h_coords;



	}

	//cout << "WTF!!1" << endl;

	is_defined = false;
	return vector3d();
}

bool quick_add_HN(const PRODART::POSE::pose_shared_ptr protein,
		const int residue_number,
		const bool overwrite){

	//cout << "called quick_add_HN for res num:" << residue_number << "\t";

	const residue_shared_ptr currRes = protein->get_residue(residue_number);
	const const_residue_shared_ptr prevRes = currRes->get_prev_residue();

	atom_shared_ptr h_0_ptr = protein->get_bb_atom(H, residue_number);

	if (h_0_ptr->isSet() && h_0_ptr->isActive() && !overwrite ){
		//cout << "HN already set\n";
		return true; //HN already set and no overwrite requested
	}

	if (currRes && prevRes && residue_number > 0){
		const_atom_shared_ptr c_m1_ptr = protein->get_bb_atom(C, residue_number - 1),
				ca_0_ptr = protein->get_bb_atom(CA, residue_number),
				n_0_ptr = protein->get_bb_atom(N, residue_number);
		if (c_m1_ptr->isSet() && ca_0_ptr->isSet() && n_0_ptr->isSet()){
			const vector3d new_h_coords = PRODART::UTILS::dihedralEnd(c_m1_ptr->get_coords(),
					ca_0_ptr->get_coords(),
					n_0_ptr->get_coords(),
					0.9992, 								// bond_3_4
					115.4200 * (PRODART::UTILS::PI / 180.0),				// angle_2_3_4
					180.0000 * (PRODART::UTILS::PI / 180.0));				// dihedral_1_2_3_4
			h_0_ptr->set_coords(new_h_coords);
			h_0_ptr->set_type(atom_type("H"));
			h_0_ptr->setSet(true);
			h_0_ptr->setActive(true);
			//cout << "!!!" << h_0_ptr->get_coords() << "!!!" << endl;
			return true;
		}
		else {
			// not enough atoms - figure out something else here

			//cout << "not enough atoms\n";


		}
	}
	else {
		// N-terminal residue - figure out something else here
		// undefined
		const_atom_shared_ptr c_0_ptr = protein->get_bb_atom(C, residue_number),
				ca_0_ptr = protein->get_bb_atom(CA, residue_number),
				n_0_ptr = protein->get_bb_atom(N, residue_number);

		if ( !c_0_ptr->isSet() || !ca_0_ptr->isSet()
				|| !n_0_ptr->isSet() ) { return false; }

		double angle = 0;

		//srand(10);
		angle = ((180.0/180.0) * PRODART::UTILS::PI);


		const vector3d new_h_coords = dihedralEnd( c_0_ptr->get_coords(),
				ca_0_ptr->get_coords(),
				n_0_ptr->get_coords(),
				0.9992,
				((115.0/180.0) * PRODART::UTILS::PI) ,
				angle );
		h_0_ptr->set_coords(new_h_coords);
		h_0_ptr->set_type(atom_type("H"));
		h_0_ptr->setSet(true);
		h_0_ptr->setActive(true);
		return true;



	}

	//cout << "WTF!!1" << endl;


	return false;
}

bool quick_add_HN(const PRODART::POSE::pose_shared_ptr protein,
		const bool overwrite){

	//cout << "started quick_add_HN to all" << endl;

	bool allOK = true;

	const int res_count = protein->get_residue_count();

	//cout << res_count << endl;

	for (int i = 0; i < res_count; i++){
		//cout << i << endl;
		bool this_result = quick_add_HN(protein, i, overwrite);
		allOK = allOK && this_result;
	}

	protein->index();
	return allOK;
}

//! quick and dirty addition of CB to a residue
bool quick_add_CB(const PRODART::POSE::pose_shared_ptr protein,
		const int residue_number,
		const bool overwrite ){

	const residue_shared_ptr currRes = protein->get_residue(residue_number);

	atom_shared_ptr CB_0_ptr = protein->get_bb_atom(CB, residue_number);

	if (CB_0_ptr->isSet()  && !overwrite ){
		//cout << "HN already set\n";
		return true; //HN already set and no overwrite requested
	}
	else {
		vector3d cb_vec = get_rough_ideal_CB(protein, residue_number);
		CB_0_ptr->set_coords(cb_vec);
		//added these two lines 28/5/12 to fix bug in pd2_motif_anneal:
		CB_0_ptr->set_type(atom_type("CB"));
		CB_0_ptr->setSet(true);
		CB_0_ptr->setActive(true);
		return true;
	}




	return false;
}

//! quick and dirty addition of  CB to a whole pose including GLYs
bool quick_add_CB_notGLY(const PRODART::POSE::pose_shared_ptr protein,
		const bool overwrite ){
	bool allOK = true;

	const int res_count = protein->get_residue_count();

	//cout << res_count << endl;

	for (int i = 0; i < res_count; i++){
		//cout << i << endl;
		if (protein->get_residue(i)->get_type().get_label3().compare("GLY") != 0){
			bool this_result = quick_add_CB(protein, i, overwrite);
			allOK = allOK && this_result;
		}
	}

	protein->index();
	return allOK;
}

bool quick_add_CB(const PRODART::POSE::pose_shared_ptr protein,
		const bool overwrite ){
	//cout << "started quick_add_HN to all" << endl;

	bool allOK = true;

	const int res_count = protein->get_residue_count();

	//cout << res_count << endl;

	for (int i = 0; i < res_count; i++){
		//cout << i << endl;
		bool this_result = quick_add_CB(protein, i, overwrite);
		allOK = allOK && this_result;
	}

	protein->index();
	return allOK;
}



UTILS::vector3d get_rough_ideal_CB(const PRODART::POSE::const_pose_shared_ptr protein,
		const int residue_number){
	if (residue_number >= protein->get_residue_count()){
		return vector3d();
	}
	if ( !protein->get_bb_atom(POSE::CA, residue_number)->isActiveAndSet()
			|| !protein->get_bb_atom(POSE::N, residue_number)->isActiveAndSet()
			|| !protein->get_bb_atom(POSE::C, residue_number)->isActiveAndSet()   ) {
		return vector3d();
	}

	const double ideal_single_residue_N_C_CA_CB_dihedral = degrees_to_radians(123.2300);	//122.6 * (UTILS::PI / 180.0); //! for non-PRO residues
	const double ideal_C_CA_CB = degrees_to_radians(111.0900);	//109.90 * (PI / 180.0); //! for non-PRO residues

	const double ideal_PRO_single_residue_N_C_CA_CB_dihedral = degrees_to_radians(113.7400);	//115.6 * (PI / 180.0);
	const double ideal_PRO_C_CA_CB = degrees_to_radians(111.7400);	//111.50 * (PI / 180.0);

	const double ideal_CA_CB_bond = 1.5461; //1.526;
	const double ideal_PRO_CA_CB_bond = 1.5399;//1.526;


	const vector3d ca = protein->get_bb_atom(POSE::CA, residue_number)->get_coords();

	const vector3d n = protein->get_bb_atom(POSE::N, residue_number)->get_coords();

	const vector3d c = protein->get_bb_atom(POSE::C, residue_number)->get_coords();

	const vector3d cb = protein->get_residue(residue_number)->get_type() == residue_type("PRO")
					?  UTILS::dihedralEnd( n, c, ca,
							ideal_PRO_CA_CB_bond,
							ideal_PRO_C_CA_CB,
							ideal_PRO_single_residue_N_C_CA_CB_dihedral )
					: UTILS::dihedralEnd( n, c, ca,
							ideal_CA_CB_bond,
							ideal_C_CA_CB,
							ideal_single_residue_N_C_CA_CB_dihedral );

	return cb;

	/* OLD method:

double ang, sx, sy;
myVector nca, cca, x, y;

myVector ca, n, c;
    double cabond = 1.5;

    Atom* cb = getAtom(CB);

ca = getAtom(CA)->getVec();

n = getAtom(N)->getVec();

c = getAtom(C)->getVec();


//vsub(ca,n,&nca);
nca = ca - n;
//vsub(ca,c,&cca);
cca = ca - c;

//vadd(nca,cca,&x);
x = nca + cca;


//vprod(nca,cca,&y);
y = nca * cca;

ang = acos(-1.0)/2.0 - asin(1.0/sqrt(3.0));
//ang = 0.9128;

sx = cabond*cos(ang)/x.mod();
sy = cabond*sin(ang)/y.mod();


cb->setxcoor( ca.x + x.x*sx + y.x*sy );
cb->setycoor( ca.y + x.y*sx + y.y*sy );
cb->setzcoor( ca.z + x.z*sx + y.z*sy );

	 */


}

void get_ca_bounding_box(const PRODART::POSE::const_pose_shared_ptr protein,
		PRODART::UTILS::vector3d& lowest, PRODART::UTILS::vector3d& highest){

	const int resCount = protein->get_residue_count();

	double x_low = std::numeric_limits<double>::max(),
		y_low = std::numeric_limits<double>::max(),
		z_low = std::numeric_limits<double>::max();

	double x_high = -std::numeric_limits<double>::max(),
		y_high = -std::numeric_limits<double>::max(),
		z_high = -std::numeric_limits<double>::max();

	for (int i = 0; i< resCount; i++){
		const_residue_shared_ptr currRes = protein->get_residue(i);

		if (currRes != 0){
			//const int numAtoms = PRODART::PROT::num_mainchain_atoms_per_residue;

				const_atom_shared_ptr currAtom = protein->get_bb_atom(CA, i);

				if (currAtom){
					if (currAtom->isActiveAndSet()){
						const vector3d coor = currAtom->get_coords();


						if (coor[0] < x_low) x_low = coor[0];
						if (coor[0] > x_high) x_high = coor[0];

						if (coor[1] < y_low) y_low = coor[1];
						if (coor[1] > y_high) y_high = coor[1];

						if (coor[2] < z_low) z_low = coor[2];
						if (coor[2] > z_high) z_high = coor[2];


					}
				}


		}

	}

	lowest[0] = x_low;
	lowest[1] = y_low;
	lowest[2] = z_low;

	highest[0] = x_high;
	highest[1] = y_high;
	highest[2] = z_high;

}

bool validate_bb_residue(const const_residue_shared_ptr res){
	bool isOK = true;
	if (res){
		residue_type ty = res->get_type();
		const bool isGLY = ty.is_equal3(residue_type("GLY"));
		if (!res->get_bb_atom(POSE::N)){
			cerr << "validate_bb_residue: ERROR: missing atom: N" << endl;
			isOK = false;
		}
		if (!res->get_bb_atom(POSE::CA)){
			cerr << "validate_bb_residue: ERROR: missing atom: CA" << endl;
			isOK = false;
		}
		if (!res->get_bb_atom(POSE::C)){
			cerr << "validate_bb_residue: ERROR: missing atom: C" << endl;
			isOK = false;
		}
		if (!res->get_bb_atom(POSE::O)){
			cerr << "validate_bb_residue: ERROR: missing atom: O" << endl;
			isOK = false;
		}
		if (!isGLY){
			if (!res->get_bb_atom(POSE::CB)){
				cerr << "validate_bb_residue: ERROR: missing atom: CB" << endl;
				isOK = false;
			}
		}
		if (isOK){
			if (!res->get_bb_atom(POSE::N)->isSet()) {
				cerr << "validate_bb_residue: ERROR: atom not set: N" << endl;
				isOK = false;
			}
			if (!res->get_bb_atom(POSE::CA)->isSet()) {
				cerr << "validate_bb_residue: ERROR: atom not set: CA" << endl;
				isOK = false;
			}
			if (!res->get_bb_atom(POSE::C)->isSet()) {
				cerr << "validate_bb_residue: ERROR: atom not set: C" << endl;
				isOK = false;
			}
			if (!res->get_bb_atom(POSE::O)->isSet()){
				cerr << "validate_bb_residue: ERROR: atom not set: O" << endl;
				isOK = false;
			}
			if (!isGLY){
				if (!res->get_bb_atom(POSE::CB)->isSet()){
					cerr << "validate_bb_residue: ERROR: atom not set: CB" << endl;
					isOK = false;
				}
			}
		}
	}
	else {
		cerr << "validate_bb_residue: ERROR: missing residue" << endl;
		isOK = false;
	}
	return isOK;
}

bool validate_bb_pose(const PRODART::POSE::const_pose_shared_ptr protein,
		const double bond_tol){
	bool isOK = true;
	protein->index();
	const int res_count = protein->get_residue_count();
	const_residue_shared_ptr currRes = const_residue_shared_ptr();
	const_chain_shared_ptr currChain = const_chain_shared_ptr();
	const_residue_shared_ptr prevRes = const_residue_shared_ptr();
	const_chain_shared_ptr prevChain = const_chain_shared_ptr();
	bool prevResOK = false;
	bool resOK = false;

	//const double bond_tol = 0.1; // was 0.09

	for (int i = 0; i < res_count; i++){
		prevResOK = resOK;
		prevRes = currRes;
		prevChain = currChain;
		currRes = protein->get_residue(i);
		currChain = currRes->get_chain();


		if (currChain->isPeptide()){
			resOK = validate_bb_residue(currRes);
			if (resOK){
				vector3d CA_0 = currRes->get_bb_atom(POSE::CA)->get_coords();
				vector3d C_0 = currRes->get_bb_atom(POSE::C)->get_coords();
				vector3d N_0 = currRes->get_bb_atom(POSE::N)->get_coords();
				if (fabs((CA_0 - C_0).mod() - 1.51) > bond_tol ){
					cerr << "validate_bb_pose: WARNING: CA_0 - C_0 dist test failed at residue\t"
							<< currRes->get_trimmed_pdb_residue_index() << "\t"
							<< currChain->getChainID() << "\t"
							<< (CA_0 - C_0).mod() << "\t"
							<< endl;
					isOK = false;
				}
				if (fabs((N_0 - CA_0).mod() - 1.46) > bond_tol ){
					cerr << "validate_bb_pose: WARNING: N_0 - CA_0 dist test failed at residue\t"
							<< currRes->get_trimmed_pdb_residue_index() << "\t"
							<< currChain->getChainID() << "\t"
							<< (N_0 - CA_0).mod() << "\t"
							<< endl;
					isOK = false;
				}

				if (prevRes && prevResOK){
					if ( prevChain == currChain){
						vector3d C_m1 = prevRes->get_bb_atom(POSE::C)->get_coords();
						if (fabs((N_0 - C_m1).mod() - 1.33) > bond_tol ){

							cerr << "validate_bb_pose: WARNING: N_0 - C_m1 dist test failed at residue\t"
									<< currRes->get_trimmed_pdb_residue_index() << "\t"
									<< currChain->getChainID() << "\t"
									<< (N_0 - C_m1).mod()
									<< endl;

							isOK = false;
						}
					}
				}
			}
			else {
				isOK = false;
				if (currRes){
					cerr << "validate_bb_pose: WARNING: at residue\t"
							<< currRes->get_trimmed_pdb_residue_index() << "\t"
							<< currChain->getChainID() << "\t"
							<< endl;
				}
				else {
					cerr << "validate_bb_pose: ERROR: residue not found\t" << endl;
				}
			}


		}
		else {
			if (currChain != prevChain){
				cerr << "validate_bb_pose: INFO: chain not peptide:\t" << currChain->getChainID() << endl;
			}
			resOK = false;
		}

	}

	return isOK;
}

bool validate_ca_pose(const PRODART::POSE::const_pose_shared_ptr protein){
	bool isOK = true;
	protein->index();
	const int res_count = protein->get_residue_count();
	const_residue_shared_ptr currRes = const_residue_shared_ptr();
	const_chain_shared_ptr currChain = const_chain_shared_ptr();
	const_residue_shared_ptr prevRes = const_residue_shared_ptr();
	const_chain_shared_ptr prevChain = const_chain_shared_ptr();
	bool prevResOK = false;
	bool resOK = false;

	for (int i = 0; i < res_count; i++){
		prevResOK = resOK;
		prevRes = currRes;
		prevChain = currChain;
		currRes = protein->get_residue(i);
		currChain = currRes->get_chain();
		if (currChain->isPeptide()){
			resOK = true;
			if (currRes){
				if (currRes->get_bb_atom(POSE::CA)){
					if (!currRes->get_bb_atom(POSE::CA)->isSet()){
						cerr << "validate_ca_pose: ERROR: atom not set: CA, at residue\t"
								<< currRes->get_trimmed_pdb_residue_index() << "\t"
								<< currChain->getChainID() << "\t"
								<< endl;
						resOK = false;
						isOK = false;
					}
				}
				else {
					cerr << "validate_ca_pose: ERROR: missing atom: CA, at residue\t"
							<< currRes->get_trimmed_pdb_residue_index() << "\t"
							<< currChain->getChainID() << "\t"
							<< endl;
					resOK = false;
					isOK = false;
				}
			}
			else {
				cerr << "validate_ca_pose: ERROR: missing residue\t"
						<< endl;
				resOK = false;
				isOK = false;
			}

			if (resOK && prevResOK){
				if ( currChain == prevChain){
					vector3d CA_0 = currRes->get_bb_atom(POSE::CA)->get_coords();
					vector3d CA_m1 = prevRes->get_bb_atom(POSE::CA)->get_coords();
					if (fabs((CA_0 - CA_m1).mod() - 3.8) > 0.3
							&& fabs((CA_0 - CA_m1).mod() - 3.0) > 0.3){
						cerr << "validate_ca_pose: ERROR: CA_0 - CA_m1 dist test failed at residue\t"
								<< currRes->get_trimmed_pdb_residue_index() << "\t"
								<< currChain->getChainID() << "\t"
								<< (CA_0 - CA_m1).mod()
								<< endl;
						isOK = false;
					}
				}

			}


		}
		else {
			if (currChain != prevChain){
				cerr << "validate_ca_pose: INFO: chain not peptide:\t" << currChain->getChainID() << endl;
			}
			resOK = false;
		}

	}
	return isOK;
}


POSE::four_state_sec_struct_vector get_phi_psi_omega_sector(POSE::const_pose_shared_ptr pose_){
	POSE::four_state_sec_struct_vector rtn_vec(pose_->get_residue_count(), POSE::ss4_UNDEF);

	for (int i = 0; i < pose_->get_residue_count(); i++){
		const double phi = pose_->get_phi(i);
		const double psi = pose_->get_psi(i);
		const double omega = pose_->get_omega_to_prev(i);
		POSE::four_state_sec_struct secs = get_phi_psi_omega_sector(phi, psi, omega);
		rtn_vec[i] = secs;
	}

	return rtn_vec;
}

void to_homopolymer(const POSE::pose_shared_ptr pose_, const POSE::residue_type rt){

	for (int i = 0 ; i < pose_->get_residue_count(); i++){
		pose_->get_residue(i)->set_type(rt);
	}
	pose_->index();

}


PRODART::POSE::residue_type_vector read_fasta_single_chain(const std::string filename){
	ifstream in_seq(filename.c_str(), std::ios::in);
	string lineStr;
	long length, lineNum = 0 ;
	string_vector SplitVec;

	PRODART::POSE::residue_type_vector rtn_seq;

	while ( !in_seq.eof() ) {
		getline(in_seq, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			if ( SplitVec[0].substr(0,1).compare(">") != 0
	                 && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() > 0 ){
				PRODART::POSE::residue_type_vector lineseq = POSE::one_letter_seq_to_residue_type_vector(lineStr);
				rtn_seq.insert(rtn_seq.end(), lineseq.begin(), lineseq.end());
			}
		}
	}

	return rtn_seq;
}

bool strip_bb_hydrogens(const POSE::pose_shared_ptr pose_){
	const std::vector<BBAtomType> bb_vec = POSE::prot_backbone_map::Instance()->get_bb_hydrogen_atom_list();
	for (int i = 0 ; i < pose_->get_residue_count(); i++){
		for (std::vector<BBAtomType>::const_iterator it = bb_vec.begin(); it != bb_vec.end(); it++){
			pose_->get_residue(i)->get_bb_atom(*it)->setSet(false);
			pose_->get_residue(i)->get_bb_atom(*it)->setActive(false);
		}
	}
	return true;
}


bool search_replace_atom_types(const POSE::pose_shared_ptr pose_, const std::map<POSE::atom_type, POSE::atom_type> s_r_map){

	/*
	std::map<POSE::atom_type, POSE::atom_type>::const_iterator iter;
	for (iter = s_r_map.begin(); iter != s_r_map.end(); iter++){
	}
	*/

	bool rtn_ok = true;

	const prot_backbone_map* bb_map = prot_backbone_map::Instance();

	for ( int index = 0; index < pose_->get_all_atom_count(); index++){
		POSE::atom_shared_ptr atm = pose_->get_atom(index);
		POSE::residue_shared_ptr res = atm->get_residue();
		POSE::atom_type atm_type = atm->get_type();
		const int rel_pos = bb_map->get_relative_location(atm_type);
		const bool is_bb  = rel_pos < bb_map->get_num_protein_backbone_atoms();
		if (s_r_map.find(atm->get_type()) != s_r_map.end()){
			POSE::atom_type rep_atm_type = s_r_map.find(atm->get_type())->second;
			const int rep_rel_pos = bb_map->get_relative_location(rep_atm_type);
			const bool rep_is_bb  = rep_rel_pos < bb_map->get_num_protein_backbone_atoms();
			if (is_bb || rep_is_bb){
				cerr << "pose_utils::search_replace_atom_types: WARNING: " << res->get_pdb_residue_index() << " " << res->get_chain()->getChainID() << " "
						<< atm_type.get_trimmed_label() << "->" << rep_atm_type.get_trimmed_label() << endl;
				cerr << "pose_utils::search_replace_atom_types: WARNING: searching and replacing with backbone atom types - results may be unpredictable at the moment" << endl;
				rtn_ok = false;
			}

			if (res->get_atom(rep_atm_type)){
				cerr << "pose_utils::search_replace_atom_types: WARNING: " << res->get_pdb_residue_index() << " " << res->get_chain()->getChainID() << " "
										<< atm_type.get_trimmed_label() << "->" << rep_atm_type.get_trimmed_label() << endl;
				cerr << "pose_utils::search_replace_atom_types: WARNING: searching and replacing will produce a duplicate atom_type" << endl;
				rtn_ok = false;
			}

			atm->set_type(rep_atm_type);

		}
	}

	return rtn_ok;
}



}
}



