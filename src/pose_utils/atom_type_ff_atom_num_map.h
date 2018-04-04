/*
 * atom_type_ff_atom_num_map.h
 *
 *  Created on: 7 Dec 2010
 *      Author: jmacdona
 */

#ifndef ATOM_TYPE_FF_ATOM_NUM_MAP_H_
#define ATOM_TYPE_FF_ATOM_NUM_MAP_H_

#include "pose/pose.h"
#include <map>

typedef std::map<PRODART::POSE::atom_type, int> atom_type_int_map;
typedef std::map<int, PRODART::POSE::atom_type> int_atom_type_map;

namespace PRODART {
namespace POSE_UTILS {

class atom_type_ff_atom_num_map {


private:



protected:


	int num_types;
	atom_type_int_map at_ff_map;
	int_atom_type_map ff_at_map;

	void add_pair(POSE::atom_type at, const int val){
		at_ff_map[at] = val;
		ff_at_map[val] = at;
		num_types++;
	}



public:

	atom_type_ff_atom_num_map() : num_types(0){
		at_ff_map.clear();
		ff_at_map.clear();
	}

	static const int NO_MATCH = 9999999;

	int get_ff_num(POSE::atom_type at) const{
		return at_ff_map.find(at) != at_ff_map.end() ? at_ff_map.find(at)->second : atom_type_ff_atom_num_map::NO_MATCH;
	}

	POSE::atom_type get_atom_type(const int val) const{
		return ff_at_map.find(val) != ff_at_map.end() ? ff_at_map.find(val)->second : POSE::atom_type();
	}



};






}
}

#endif /* ATOM_TYPE_FF_ATOM_NUM_MAP_H_ */
