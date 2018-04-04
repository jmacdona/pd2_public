//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * sim_cell.h
 *
 *  Created on: 4 Jun 2010
 *      Author: jmacdona
 */

#ifndef SIM_CELL_H_
#define SIM_CELL_H_


#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <map>
#include <list>
#include <iomanip>

#include "pose/pose.h"
#include "pose_meta_defs.h"



namespace PRODART {
namespace POSE {
namespace META {




class sim_cell {

public:

	sim_cell();
	~sim_cell();

	void set_cutoff(double dist);

	void set_x_limits(double min, double max);
	void set_y_limits(double min, double max);
	void set_z_limits(double min, double max);
	//void set_num_residues(const int num);
	void set_num_atoms(const int num);

	//! for debugging
	void print_vectors();

	int get_cell_num(int xgrid, int ygrid, int zgrid);


	void clear_cells(void);
	void clear_atom_list();

	bool add_to_cell(const PRODART::POSE::const_atom_shared_ptr atm);

	//! for debugging only
	/*
	void print_cell_freqs(){
		sim_sub_cell_vector::iterator iter;
		int cell_num = 0;
		int sum = 0;
		for (iter = cells.begin(); iter != cells.end(); iter++){
			cout << cell_num << "\t" << iter->atom_list.size() << endl;
			sum += iter->atom_list.size();
			cell_num++;
		}
		cout << "sum:\t" << sum << endl;
	}
	*/

	bool create_cells();
	bool create_cell_neighbour_vec();


	void recalcPairListDists(nb_ele_vector& nb_pair_list);
	void updatePairList(nb_ele_vector& nb_pair_list);


	//! trims pairs below cutoff
	void updatePairList_trim(nb_ele_vector& nb_pair_list,
			const_pose_shared_ptr protein);



protected:




	inline void addInterCellPairsToList(nb_ele_vector& nb_pair_list, const int num_cell1, const int num_cell2);
	inline void addIntraCellPairsToList(nb_ele_vector& nb_pair_list, const int num_cell);



	inline void addInterCellPairsToList_trim(nb_ele_vector& nb_pair_list,
			const int num_cell1,
			const int num_cell2,
			const_pose_shared_ptr protein);
	inline void addIntraCellPairsToList_trim(nb_ele_vector& nb_pair_list,
			const int num_cell,
			const_pose_shared_ptr protein);




	double cutoff;

	double x_min, x_max,
		y_min, y_max,
		z_min, z_max;

	int x_max_grid,
		y_max_grid,
		z_max_grid;

	int num_atoms;

	int_pair_vector cell_neighbour_vec;
	int_vector head;
	int_vector lscl;

	PRODART::POSE::const_atom_shared_ptr_vector atom_list;

	static const int LEMPTY;

};

/* *************************************
 *  INLINE functions
 * *************************************/




inline int sim_cell::get_cell_num(int xgrid, int ygrid, int zgrid){
	return xgrid + (x_max_grid * ygrid) + (x_max_grid * y_max_grid * zgrid);
}





inline void sim_cell::clear_cells(void){
	/*
	sim_sub_cell_vector::iterator iter;
	for (iter = cells.begin(); iter != cells.end(); iter++){
		iter->atom_list.clear();
	}
	*/

	int_vector::iterator int_iter;

	for (int_iter = head.begin(); int_iter != head.end(); int_iter++){
		*int_iter = LEMPTY;
	}

	for (int_iter = lscl.begin(); int_iter != lscl.end(); int_iter++){
		*int_iter = LEMPTY;
	}

	this->clear_atom_list();

}

inline void sim_cell::clear_atom_list(){
	//this->atom_list.clear();
	this->atom_list.resize(0);
	if (int(atom_list.capacity()) < num_atoms){
		this->atom_list.reserve(num_atoms);
	}
}


inline bool sim_cell::add_to_cell(const PRODART::POSE::const_atom_shared_ptr atm){

	this->atom_list.push_back(atm);
	const PRODART::UTILS::vector3d pos = atm->get_coords();
	const int atom_index = atom_list.size() - 1;
	const double x_val = pos.x - x_min;
	const int xgrid = static_cast<int>((x_val) / cutoff);
	if (x_val < 0 || xgrid < 0 || xgrid >= x_max_grid) return false;
	const double y_val = pos.y - y_min;
	const int ygrid = static_cast<int>((y_val) / cutoff);
	if (y_val < 0 || ygrid < 0 || ygrid >= y_max_grid) return false;
	const double z_val = pos.z - z_min;
	const int zgrid = static_cast<int>((z_val) / cutoff);
	if (z_val < 0 || zgrid < 0 || zgrid >= z_max_grid) return false;

	const int cell_num = this->get_cell_num(xgrid, ygrid, zgrid);

	//link to previous head
	lscl[atom_index] = head[cell_num];
	//current becomes new head
	head[cell_num] = atom_index;

	//cells[cell_num].atom_list.push_back(atom_index);

	return true;

}



inline bool sim_cell::create_cells(){

	x_max_grid = static_cast<int>((x_max - x_min) / cutoff) + 1;
	y_max_grid = static_cast<int>((y_max - y_min) / cutoff) + 1;
	z_max_grid = static_cast<int>((z_max - z_min) / cutoff) + 1;

	const int num_cells = x_max_grid * y_max_grid * z_max_grid;

	//head.clear();
	//lscl.clear();
	head.resize(0);
	lscl.resize(0);
	head.resize(num_cells, LEMPTY);
	lscl.resize(num_atoms + 100, LEMPTY);
	atom_list.clear();
	atom_list.reserve(num_atoms+100);

	//cells.resize(num_cells);

	return true;
}




}
}
}



/*
 *
 *
 *
 *

	inline void addInterCellPairsToList_for_others(nb_ele_vector& nb_pair_list, const int num_cell1, const int num_cell2);
	inline void addIntraCellPairsToList_for_others(nb_ele_vector& nb_pair_list, const int num_cell);

	inline void addInterCellPairsToList_for_others_trim(nb_ele_vector& nb_pair_list,
			const int num_cell1,
			const int num_cell2,
			const_pose_shared_ptr protein);
	inline void addIntraCellPairsToList_for_others_trim(nb_ele_vector& nb_pair_list,
			const int num_cell,
			const_pose_shared_ptr protein);


	//! trims pairs below cutoff
	void updatePairList_for_others_trim(nb_ele_vector& nb_pair_list,
			const const_pose_shared_ptr protein);


inline void addInterCellPairsToList(nb_ele_list& nb_pair_list, const int num_cell1, const int num_cell2);
inline void addIntraCellPairsToList(nb_ele_list& nb_pair_list, const int num_cell);

inline void addInterCellPairsToList_for_others(nb_ele_list& nb_pair_list, const int num_cell1, const int num_cell2);
inline void addIntraCellPairsToList_for_others(nb_ele_list& nb_pair_list, const int num_cell);
 */


#endif /* SIM_CELL_H_ */
