/*
 * atom_type_map.h
 *
 *  Created on: 1 Aug 2011
 *      Author: jmacdona
 */

#ifndef ATOM_TYPE_MAP_H_
#define ATOM_TYPE_MAP_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>
#include <vector>
#include <string>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/thread.hpp>

#include "prodart_env/prodart_env.h"

#include "atom_type.h"

namespace PRODART {
namespace POSE {

typedef std::map<std::string, int> string_int_map;
typedef std::map<int, std::string> int_string_map;

//! singleton class
class atom_type_map{

private:
	atom_type_map(const atom_type_map&);


	static void Init();

	string_int_map str_to_ff_id_map;
	int_string_map ff_id_to_str_map;

protected:
	atom_type_map();
	virtual ~atom_type_map();

public:

	static const atom_type_map* Instance();


	int get_ff_id(const atom_type &at) const;


};


inline int atom_type_map::get_ff_id(const atom_type &at) const{
	return str_to_ff_id_map.find(at.get_label()) != str_to_ff_id_map.end() ? str_to_ff_id_map.find(at.get_label())->second : atom_type::UNK_TYPE;
}



}
}

#endif /* ATOM_TYPE_MAP_H_ */
