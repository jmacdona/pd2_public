/*
 * atom_type.cpp
 *
 *  Created on: 1 Aug 2011
 *      Author: jmacdona
 */

#include "atom_type.h"
#include "atom_type_map.h"

using boost::trim;

namespace PRODART {
namespace POSE {


void atom_type::auto_set_ff_atom_id(){
	const atom_type_map* at_map = atom_type_map::Instance();
	this->set_ff_atom_id(at_map->get_ff_id(*this));
}

std::string atom_type::get_trimmed_label() const{
	std::string ret_str = this->get_label();
	trim(ret_str);
	return ret_str;
}


}
}
