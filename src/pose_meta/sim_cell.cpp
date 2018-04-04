//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * sim_cell.cpp
 *
 *  Created on: 4 Jun 2010
 *      Author: jmacdona
 */
#include "sim_cell.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::UTILS;


namespace PRODART {
namespace POSE {
namespace META {







const int sim_cell::LEMPTY = -1;

sim_cell::sim_cell(){

	cutoff = 15;

	x_min = -200;
	x_max = 200;
	y_min = -200;
	y_max = 200;
	z_min = -200;
	z_max = 200;

	num_atoms = 0;

	x_max_grid = static_cast<int>((x_min - x_max) / cutoff) + 1;
	y_max_grid = static_cast<int>((y_min - y_max) / cutoff) + 1;
	z_max_grid = static_cast<int>((z_min - z_max) / cutoff) + 1;

	//cells.clear();

	head.clear();
	lscl.clear();
	atom_list.clear();

}

sim_cell::~sim_cell(){

}


void sim_cell::set_cutoff(double dist){
	cutoff = dist;
}

void sim_cell::set_x_limits(double min, double max){
	x_min = min;
	x_max = max;
}
void sim_cell::set_y_limits(double min, double max){
	y_min = min;
	y_max = max;
}
void sim_cell::set_z_limits(double min, double max){
	z_min = min;
	z_max = max;
}

void sim_cell::set_num_atoms(const int num){
	num_atoms = num;
	atom_list.reserve(num_atoms);
}


bool sim_cell::create_cell_neighbour_vec(){
	cell_neighbour_vec.resize(0);
	const unsigned int ressize = head.size() * 14;
	if (cell_neighbour_vec.capacity() < ressize){
		cell_neighbour_vec.reserve(ressize);
	}
	/*
	cout << "maximums:\t"
		 << x_max_grid << " "
		 << y_max_grid << " "
		 << z_max_grid << "\t"
		 << endl;
	 */

	for (int x = 0; x < x_max_grid; x++ ){
		for (int y = 0; y < y_max_grid; y++ ){
			for (int z = 0; z < z_max_grid; z++ ){



				const int this_cell = this->get_cell_num(x,y,z);

				for (int x_add = x - 1; x_add <= x + 1; x_add++ ){
					for (int y_add = y - 1; y_add <= y + 1; y_add++ ){
						for (int z_add = z - 1; z_add <= z + 1; z_add++ ){

							const int neigh_cell = this->get_cell_num(x_add, y_add, z_add);

							if (neigh_cell > this_cell
									&& x_add < x_max_grid
									&& y_add < y_max_grid
									&& z_add < z_max_grid
									&& x_add >= 0
									&& y_add >= 0
									&& z_add >= 0){
								/*
								cout << x << " "
									 << y << " "
									 << z << "\t"
									 << this_cell << "\t" << neigh_cell << endl;
								*/

								cell_neighbour_vec.push_back(int_pair(this_cell, neigh_cell));

							}

						}
					}
				}



			}
		}
	}
	return true;
}






void sim_cell::recalcPairListDists(nb_ele_vector& nb_pair_list){
	nb_ele_vector::iterator iter;

	for (iter = nb_pair_list.begin(); iter != nb_pair_list.end(); iter++){
		if (iter->atom1_ptr->has_moved() || iter->atom2_ptr->has_moved()){
			const vector3d tempVec_i = iter->atom1_ptr->get_coords(); //protein.getVecByAtomIndex(i);
			const vector3d tempVec_j = iter->atom2_ptr->get_coords(); //protein.getVecByAtomIndex(j);
			iter->dist_sq = (tempVec_i - tempVec_j).mod_sq();
			iter->dist = sqrt(iter->dist_sq);
		}
	}
}

void sim_cell::updatePairList(nb_ele_vector& nb_pair_list){
	nb_pair_list.reserve(lscl.size()*400);
	int_pair_vector::iterator iter;
	//cout << "inter pairs" << endl;
	for (iter = cell_neighbour_vec.begin(); iter != cell_neighbour_vec.end(); iter++){
		this->addInterCellPairsToList(nb_pair_list, iter->first, iter->second);
	}

	//sim_sub_cell_vector::iterator intra_it;
	//cout << "intra pairs" << endl;
	const int num_cells = static_cast<int>(head.size());
	for(int cn = 0; cn < num_cells; cn++){
		this->addIntraCellPairsToList(nb_pair_list, cn);
	}


}




void sim_cell::updatePairList_trim(nb_ele_vector& nb_pair_list,
		const_pose_shared_ptr protein){
	//cout << "reserving\t" << lscl.size()*400 << endl;
	const unsigned int ressize =  lscl.size()*400;
	if (nb_pair_list.capacity() < ressize){
		nb_pair_list.reserve(ressize);
	}
	//cout << "reserved\t" << lscl.size()*400 << endl;
	int_pair_vector::iterator iter;
	//cout << "inter pairs" << endl;
	for (iter = cell_neighbour_vec.begin(); iter != cell_neighbour_vec.end(); iter++){
		this->addInterCellPairsToList_trim(nb_pair_list,
				iter->first,
				iter->second,
				protein);
	}

	//sim_sub_cell_vector::iterator intra_it;
	//cout << "intra pairs" << endl;
	const int num_cells = static_cast<int>(head.size());
	for(int cn = 0; cn < num_cells; cn++){
		this->addIntraCellPairsToList_trim(nb_pair_list,
				cn,
				protein);
	}


}





void sim_cell::print_vectors(){

	int_vector::iterator int_iter;
	cout << "head:\n";
	int i = 0;
	for (int_iter = head.begin(); int_iter != head.end(); int_iter++){
		cout << i << ":" << *int_iter << "\t";
		i++;
	}

	cout << "\nlscl:\n";
	i = 0;
	for (int_iter = lscl.begin(); int_iter != lscl.end(); int_iter++){
		cout << i << ":" <<  *int_iter << "\t";
		i++;
	}

	cout << endl;


}


inline void sim_cell::addInterCellPairsToList(nb_ele_vector& nb_pair_list, const int num_cell1, const int num_cell2){

	int atm_c1 = head[num_cell1];
	int atm_c2 = head[num_cell2];

	if (atm_c1 == LEMPTY || atm_c2 == LEMPTY) return;

	while (atm_c1 != LEMPTY){
		atm_c2 = head[num_cell2];
		while (atm_c2 != LEMPTY){

			nb_pair_element b_ele;
			b_ele.atom1_ptr = this->atom_list[atm_c1];
			b_ele.atom2_ptr = this->atom_list[atm_c2];
			b_ele.atom1 = this->atom_list[atm_c1]->get_seq_num(); //atm_c1;
			b_ele.atom2 = this->atom_list[atm_c2]->get_seq_num(); //atm_c2;
			b_ele.atype1 = this->atom_list[atm_c1]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c1);
			b_ele.atype2 = this->atom_list[atm_c2]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c2);

			const vector3d tempVec_i = this->atom_list[atm_c1]->get_coords(); //protein.getVecByAtomIndex(i);
			const vector3d tempVec_j = this->atom_list[atm_c2]->get_coords(); //protein.getVecByAtomIndex(j);
			b_ele.dist_sq = (tempVec_i - tempVec_j).mod_sq();
			b_ele.dist = sqrt(b_ele.dist_sq);
			nb_pair_list.push_back(b_ele);

			atm_c2 = lscl[atm_c2];
		}
		atm_c1 = lscl[atm_c1];
	}
}

inline void sim_cell::addIntraCellPairsToList(nb_ele_vector& nb_pair_list, const int num_cell){
	int atm_c1 = head[num_cell];
	while (atm_c1 != LEMPTY){
		int atm_c2 = lscl[atm_c1];
		while (atm_c2 != LEMPTY){

			nb_pair_element b_ele;
			b_ele.atom1_ptr = this->atom_list[atm_c1];
			b_ele.atom2_ptr = this->atom_list[atm_c2];
			b_ele.atom1 = this->atom_list[atm_c1]->get_seq_num(); //atm_c1;
			b_ele.atom2 = this->atom_list[atm_c2]->get_seq_num(); //atm_c2;
			b_ele.atype1 = this->atom_list[atm_c1]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c1);
			b_ele.atype2 = this->atom_list[atm_c2]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c2);

			const vector3d tempVec_i = this->atom_list[atm_c1]->get_coords(); //protein.getVecByAtomIndex(i);
			const vector3d tempVec_j = this->atom_list[atm_c2]->get_coords(); //protein.getVecByAtomIndex(j);
			b_ele.dist_sq = (tempVec_i - tempVec_j).mod_sq();
			b_ele.dist = sqrt(b_ele.dist_sq);
			nb_pair_list.push_back(b_ele);

			atm_c2 = lscl[atm_c2];
		}
		atm_c1 = lscl[atm_c1];
	}


}







inline void sim_cell::addInterCellPairsToList_trim(nb_ele_vector& nb_pair_list,
		const int num_cell1,
		const int num_cell2,
		const_pose_shared_ptr protein){

	int atm_c1 = head[num_cell1];
	int atm_c2 = head[num_cell2];

	if (atm_c1 == LEMPTY || atm_c2 == LEMPTY) return;

	while (atm_c1 != LEMPTY){
		atm_c2 = head[num_cell2];
		while (atm_c2 != LEMPTY){



			//const int i = b_ele.atom1;
			//const int j = b_ele.atom2;
			const vector3d tempVec_i = this->atom_list[atm_c1]->get_coords(); //protein.getVecByAtomIndex(i);
			const vector3d tempVec_j = this->atom_list[atm_c2]->get_coords(); //protein.getVecByAtomIndex(j);
			const double dist_sq = (tempVec_i - tempVec_j).mod_sq();
			const double dist = sqrt(dist_sq);
			if (dist < this->cutoff){
				nb_pair_element b_ele;
				b_ele.atom1_ptr = this->atom_list[atm_c1];
				b_ele.atom2_ptr = this->atom_list[atm_c2];
				b_ele.atom1 = this->atom_list[atm_c1]->get_seq_num(); //atm_c1;
				b_ele.atom2 = this->atom_list[atm_c2]->get_seq_num(); //atm_c2;
				b_ele.atype1 = this->atom_list[atm_c1]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c1);
				b_ele.atype2 = this->atom_list[atm_c2]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c2);
				b_ele.dist_sq = dist_sq;
				b_ele.dist = dist;
				b_ele.res_num1 = this->atom_list[atm_c1]->get_residue()->get_internal_residue_index();
				b_ele.res_num2 = this->atom_list[atm_c2]->get_residue()->get_internal_residue_index();

				nb_pair_list.push_back(b_ele);
			}
			atm_c2 = lscl[atm_c2];
		}
		atm_c1 = lscl[atm_c1];
	}
}

inline void sim_cell::addIntraCellPairsToList_trim(nb_ele_vector& nb_pair_list,
		const int num_cell,
		const_pose_shared_ptr protein){
	int atm_c1 = head[num_cell];
	while (atm_c1 != LEMPTY){
		int atm_c2 = lscl[atm_c1];
		while (atm_c2 != LEMPTY){



			//const int i = b_ele.atom1;
			//const int j = b_ele.atom2;
			const vector3d tempVec_i = this->atom_list[atm_c1]->get_coords(); //protein.getVecByAtomIndex(i);
			const vector3d tempVec_j = this->atom_list[atm_c2]->get_coords(); //protein.getVecByAtomIndex(j);
			const double dist_sq = (tempVec_i - tempVec_j).mod_sq();
			const double dist = sqrt(dist_sq);
			if (dist < this->cutoff){
				nb_pair_element b_ele;
				b_ele.atom1_ptr = this->atom_list[atm_c1];
				b_ele.atom2_ptr = this->atom_list[atm_c2];
				b_ele.atom1 = this->atom_list[atm_c1]->get_seq_num(); //atm_c1;
				b_ele.atom2 = this->atom_list[atm_c2]->get_seq_num(); //atm_c2;
				b_ele.atype1 = this->atom_list[atm_c1]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c1);
				b_ele.atype2 = this->atom_list[atm_c2]->get_type(); //PRODART::PROT::getAtomTypeFromIndex(atm_c2);
				b_ele.dist_sq = dist_sq;
				b_ele.dist = dist;
				b_ele.res_num1 = this->atom_list[atm_c1]->get_residue()->get_internal_residue_index();
				b_ele.res_num2 = this->atom_list[atm_c2]->get_residue()->get_internal_residue_index();

				nb_pair_list.push_back(b_ele);
			}

			atm_c2 = lscl[atm_c2];
		}
		atm_c1 = lscl[atm_c1];
	}


}
























}
}
}









/*
 * JUNK


inline void sim_cell::addInterCellPairsToList_for_others_trim(nb_ele_vector& nb_pair_list,
		const int num_cell1,
		const int num_cell2,
		const_pose_shared_ptr protein){

	int res_c1 = head[num_cell1];
	int res_c2 = head[num_cell2];

	if (res_c1 == LEMPTY || res_c2 == LEMPTY) return;

	while (res_c1 != LEMPTY){
		res_c2 = head[num_cell2];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			const PRODART::PROT::AtomType at_i = b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			const PRODART::PROT::AtomType at_j = b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);

			if (!(((at_i == PRODART::PROT::H || at_i == PRODART::PROT::O)&&(at_j == PRODART::PROT::H || at_j == PRODART::PROT::O))
					|| (at_i == PRODART::PROT::CA && at_j == PRODART::PROT::CA)
					|| (at_i == PRODART::PROT::CB && at_j == PRODART::PROT::CB))){

				const int i = b_ele.atom1;
				const int j = b_ele.atom2;
				const myVector tempVec_i = protein.getVecByAtomIndex(i);
				const myVector tempVec_j = protein.getVecByAtomIndex(j);
				b_ele.dist_sq = (tempVec_i - tempVec_j).mod_sq();
				b_ele.dist = sqrt(b_ele.dist_sq);
				if (b_ele.dist < this->cutoff){
					nb_pair_list.push_back(b_ele);
				}
			}

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}
}

inline void sim_cell::addIntraCellPairsToList_for_others_trim(nb_ele_vector& nb_pair_list,
		const int num_cell,
		const_pose_shared_ptr protein){

	int res_c1 = head[num_cell];
	while (res_c1 != LEMPTY){
		int res_c2 = lscl[res_c1];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			const PRODART::PROT::AtomType at_i = b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			const PRODART::PROT::AtomType at_j = b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);

			if (!(((at_i == PRODART::PROT::H || at_i == PRODART::PROT::O)&&(at_j == PRODART::PROT::H || at_j == PRODART::PROT::O))
					|| (at_i == PRODART::PROT::CA && at_j == PRODART::PROT::CA)
					|| (at_i == PRODART::PROT::CB && at_j == PRODART::PROT::CB))){

				const int i = b_ele.atom1;
				const int j = b_ele.atom2;
				const myVector tempVec_i = protein.getVecByAtomIndex(i);
				const myVector tempVec_j = protein.getVecByAtomIndex(j);
				b_ele.dist_sq = (tempVec_i - tempVec_j).mod_sq();
				b_ele.dist = sqrt(b_ele.dist_sq);
				if (b_ele.dist < this->cutoff){
					nb_pair_list.push_back(b_ele);
				}
			}

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}


}


inline void sim_cell::addInterCellPairsToList_for_others(nb_ele_vector& nb_pair_list, const int num_cell1, const int num_cell2){

	int res_c1 = head[num_cell1];
	int res_c2 = head[num_cell2];

	if (res_c1 == LEMPTY || res_c2 == LEMPTY) return;

	while (res_c1 != LEMPTY){
		res_c2 = head[num_cell2];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			const PRODART::PROT::AtomType at_i = b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			const PRODART::PROT::AtomType at_j = b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);

			if (!(((at_i == PRODART::PROT::H || at_i == PRODART::PROT::O)&&(at_j == PRODART::PROT::H || at_j == PRODART::PROT::O))
					|| (at_i == PRODART::PROT::CA && at_j == PRODART::PROT::CA)
					|| (at_i == PRODART::PROT::CB && at_j == PRODART::PROT::CB))){
				nb_pair_list.push_back(b_ele);
			}

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}
}

inline void sim_cell::addIntraCellPairsToList_for_others(nb_ele_vector& nb_pair_list, const int num_cell){
	int res_c1 = head[num_cell];
	while (res_c1 != LEMPTY){
		int res_c2 = lscl[res_c1];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			const PRODART::PROT::AtomType at_i = b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			const PRODART::PROT::AtomType at_j = b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);

			if (!(((at_i == PRODART::PROT::H || at_i == PRODART::PROT::O)&&(at_j == PRODART::PROT::H || at_j == PRODART::PROT::O))
					|| (at_i == PRODART::PROT::CA && at_j == PRODART::PROT::CA)
					|| (at_i == PRODART::PROT::CB && at_j == PRODART::PROT::CB))){
				nb_pair_list.push_back(b_ele);
			}

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}


}

inline void sim_cell::addInterCellPairsToList(nb_ele_list& nb_pair_list, const int num_cell1, const int num_cell2){

	int res_c1 = head[num_cell1];
	int res_c2 = head[num_cell2];

	if (res_c1 == LEMPTY || res_c2 == LEMPTY) return;

	while (res_c1 != LEMPTY){
		res_c2 = head[num_cell2];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);
			nb_pair_list.push_back(b_ele);

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}
}

inline void sim_cell::addIntraCellPairsToList(nb_ele_list& nb_pair_list, const int num_cell){
	int res_c1 = head[num_cell];
	while (res_c1 != LEMPTY){
		int res_c2 = lscl[res_c1];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);
			nb_pair_list.push_back(b_ele);

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}


}

inline void sim_cell::addInterCellPairsToList_for_others(nb_ele_list& nb_pair_list, const int num_cell1, const int num_cell2){

	int res_c1 = head[num_cell1];
	int res_c2 = head[num_cell2];

	if (res_c1 == LEMPTY || res_c2 == LEMPTY) return;

	while (res_c1 != LEMPTY){
		res_c2 = head[num_cell2];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			const PRODART::PROT::AtomType at_i = b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			const PRODART::PROT::AtomType at_j = b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);

			if (!(((at_i == PRODART::PROT::H || at_i == PRODART::PROT::O)&&(at_j == PRODART::PROT::H || at_j == PRODART::PROT::O))
					|| (at_i == PRODART::PROT::CA && at_j == PRODART::PROT::CA)
					|| (at_i == PRODART::PROT::CB && at_j == PRODART::PROT::CB))){
				nb_pair_list.push_back(b_ele);
			}

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}
}

inline void sim_cell::addIntraCellPairsToList_for_others(nb_ele_list& nb_pair_list, const int num_cell){
	int res_c1 = head[num_cell];
	while (res_c1 != LEMPTY){
		int res_c2 = lscl[res_c1];
		while (res_c2 != LEMPTY){

			t_nb_pair_element b_ele;
			b_ele.atom1 = res_c1;
			b_ele.atom2 = res_c2;
			const PRODART::PROT::AtomType at_i = b_ele.atype1 = PRODART::PROT::getAtomTypeFromIndex(res_c1);
			const PRODART::PROT::AtomType at_j = b_ele.atype2 = PRODART::PROT::getAtomTypeFromIndex(res_c2);

			if (!(((at_i == PRODART::PROT::H || at_i == PRODART::PROT::O)&&(at_j == PRODART::PROT::H || at_j == PRODART::PROT::O))
					|| (at_i == PRODART::PROT::CA && at_j == PRODART::PROT::CA)
					|| (at_i == PRODART::PROT::CB && at_j == PRODART::PROT::CB))){
				nb_pair_list.push_back(b_ele);
			}

			res_c2 = lscl[res_c2];
		}
		res_c1 = lscl[res_c1];
	}


}


	inline void old_addInterCellPairsToList(BondList_vector& nb_pair_list, const int num_cell1, const int num_cell2){

		sim_sub_cell& sub_cell1 = cells[num_cell1];
		sim_sub_cell& sub_cell2 = cells[num_cell2];

		if (sub_cell1.atom_list.size() == 0 || sub_cell1.atom_list.size() == 0){
			return;
		}
		std::list<int>::const_iterator iter1;
		for (iter1 = sub_cell1.atom_list.begin(); iter1 != sub_cell1.atom_list.end(); iter1++){
			std::list<int>::const_iterator iter2;
			for (iter2 = sub_cell2.atom_list.begin(); iter2 != sub_cell2.atom_list.end(); iter2++){

				t_Bond_list_element b_ele;
				b_ele.atom1 = *iter1;
				b_ele.atom2 = *iter2;
				nb_pair_list.push_back(b_ele);

			}
		}

	}

	inline void old_addIntraCellPairsToList(BondList_vector& nb_pair_list, const int num_cell){
		sim_sub_cell& sub_cell = cells[num_cell];
		if (sub_cell.atom_list.size() == 0 ){
			return;
		}
		std::list<int>::const_iterator iter1;
		std::list<int>::const_iterator iter2;
		for (iter1 = sub_cell.atom_list.begin(); iter1 != sub_cell.atom_list.end(); iter1++){
			for (iter2 = iter1; iter2 != sub_cell.atom_list.end(); iter2++){
				if(iter1 != iter2){
					t_Bond_list_element b_ele;
					b_ele.atom1 = *iter1;
					b_ele.atom2 = *iter2;
					nb_pair_list.push_back(b_ele);
				}
			}
		}

	}
 */





