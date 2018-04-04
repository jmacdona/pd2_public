//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * atom_residue_types.h
 *
 *  Created on: Jan 28, 2010
 *      Author: jmacdon
 */

#ifndef ATOM_RESIDUE_TYPES_H_
#define ATOM_RESIDUE_TYPES_H_

#include <string>
#include <vector>
#include <iostream>
#include <cassert>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/functional/hash.hpp>


namespace PRODART {
namespace POSE {

class residue_type;

typedef std::vector<residue_type> residue_type_vector;


class residue_type{

	//! comparison using atom_label[4]
	friend bool operator==( const residue_type&,  const residue_type& );
	//! comparison using atom_label[4]
	friend bool operator<( const residue_type&,  const residue_type& );

	friend std::size_t hash_value(const residue_type& type);

	static const unsigned int max_len = 5;
	static const unsigned int pdb_format_len = 3;

public:

	residue_type();
	residue_type(const std::string label);
	residue_type& operator=(const residue_type&);

	std::string get_label() const;
	std::string get_label3() const;
	std::string get_pdb_formatted_label() const;

	//compare first 3 characters of label
	bool is_equal3(const residue_type&) const;

	//std::string get_ff_string() const;

	int get_ff_res_id() const;

private:

	void set_ff_res_id(const int);

	//! name as read/output to PDB file and used to index atom in residue so must be unique.
	char residue_label[max_len];
	//! forcefield residue id
	int ff_res_id;
	//! residue type as defined in force field - may not be the same as atom_name above (CAUTION: may use different pointer based formulation)
	//char ff_atom_type[4];

};




/*
 * ************** INLINE FUNCTIONS ************
 */

inline std::size_t hash_value(const residue_type& type){

	boost::hash<char[residue_type::max_len]> hasher;
	return hasher(type.residue_label);


	//return static_cast<std::size_t>(type.get_ff_atom_id());
}


inline residue_type::residue_type(){
	for (unsigned int i = 0; i < max_len; i++){
		this->residue_label[i] = ' ';
	}
	ff_res_id = 0;
}

inline residue_type::residue_type(const std::string label){
	const unsigned int str_len = label.length();
	for (unsigned int i = 0; i < max_len; i++){
		if (i < str_len){
			this->residue_label[i] = label[i];
		}
		else {
			this->residue_label[i] = ' ';
		}
	}
	ff_res_id = 0;

	if (str_len > max_len){
		std::cerr << "residue_type: WARNING: label exceeds max size" << std::endl;
		assert(str_len > max_len);
	}

}

inline residue_type& residue_type::operator=(const residue_type& rt){
	for (unsigned int i = 0; i < max_len; i++){
		this->residue_label[i] = rt.residue_label[i];
	}
	this->ff_res_id = rt.ff_res_id;
	return *this;
}

inline std::string residue_type::get_label() const{
	std::string ret_str(max_len, ' ');
	for (unsigned int i = 0; i < max_len; i++){
		ret_str[i] = this->residue_label[i];
	}
	return ret_str;
}

inline std::string residue_type::get_label3() const{
	std::string ret_str(3, ' ');
	for (unsigned int i = 0; i < 3; i++){
		ret_str[i] = this->residue_label[i];
	}
	return ret_str;
}

inline std::string residue_type::get_pdb_formatted_label() const{
	std::string ret_str = this->get_label();
	boost::trim(ret_str);
	const unsigned int size = ret_str.size();

	if (size == pdb_format_len){
		return ret_str;
	}
	else if (size > pdb_format_len){
		ret_str = ret_str.substr(0,pdb_format_len);
	}
	else if (size < pdb_format_len) {
		//size = ret_str.size();
		for (unsigned int i = size; i < pdb_format_len; i++){
			ret_str.insert(0," ");
		}
	}
	else if (size == 0){
		return std::string(pdb_format_len, ' ');;//std::string("   ");
	}
	return ret_str;
}

inline void residue_type::set_ff_res_id(const int val){
	this->ff_res_id = val;
}

inline int residue_type::get_ff_res_id() const{
	return this->ff_res_id;
}

inline bool residue_type::is_equal3(const residue_type& rt2) const{
	for (unsigned int i = 0; i < 3; i++){
		if (this->residue_label[i] != rt2.residue_label[i]){
			return false;
		}
	}
	return true;
}

/*
inline std::string residue_type::get_ff_string() const{
	return this->get_label();
}
*/

inline bool operator==( const residue_type& rt1,  const residue_type& rt2 ){
	for (unsigned int i = 0; i < residue_type::max_len; i++){
		if (rt1.residue_label[i] != rt2.residue_label[i]){
			return false;
		}
	}
	return true;
}

inline bool operator<( const residue_type& rt1,  const residue_type& rt2 ){
	for (unsigned int i = 0; i < residue_type::max_len; i++){
		if (rt1.residue_label[i] < rt2.residue_label[i]){
			return true;
		}
		else if (rt1.residue_label[i] > rt2.residue_label[i]) {
			return false;
		}
	}
	return false;
}


}
}

#endif /* ATOM_RESIDUE_TYPES_H_ */
