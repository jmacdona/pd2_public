//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * atom_type.h
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#ifndef ATOM_TYPE_H_
#define ATOM_TYPE_H_

#include <string>
#include <cctype>
#include <vector>
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/functional/hash.hpp>

namespace PRODART {
namespace POSE {

class atom_type;



typedef boost::tuple<POSE::atom_type, POSE::atom_type> atom_type2_tuple;


class atom_type{

	//! comparison using atom_label[4]
	friend bool operator==( const atom_type&,  const atom_type& );
	//! comparison using atom_label[4]
	friend bool operator!=( const atom_type&,  const atom_type& );
	//! comparison using atom_label[4]
	friend bool operator<( const atom_type&,  const atom_type& );

	friend std::size_t hash_value(const atom_type2_tuple& tup);
	friend std::size_t hash_value(const atom_type& type);

public:

	static const int UNK_TYPE = -1;

	static const unsigned int MAX_LABEL_LEN = 4;

	atom_type();
	atom_type(const std::string label, const bool noff = false);
	atom_type& operator=(const atom_type&);


	std::string get_label() const;
	std::string get_trimmed_label() const;
	std::string get_pdb_formatted_label() const;
	//std::string get_ff_string() const;

	int get_ff_atom_id() const;
	void auto_set_ff_atom_id();

private:

	void set_ff_atom_id(const int);


	//! name as read/output to PDB file and used to index atom in residue so must be unique.
	char atom_label[MAX_LABEL_LEN];
	//! forcefield atom id
	int ff_atom_id;
	//! atom type as defined in force field - may not be the same as atom_name above (CAUTION: may use different pointer based formulation)
	//char ff_atom_type[4];

};

typedef std::vector<atom_type>  atom_type_vector;
typedef std::vector<atom_type_vector>  atom_type_vector_vector;


/*
 * ************** INLINE FUNCTIONS ************
 */


inline std::size_t hash_value(const atom_type& type){

	boost::hash<char[atom_type::MAX_LABEL_LEN]> hasher;
	return hasher(type.atom_label);


	//return static_cast<std::size_t>(type.get_ff_atom_id());
}

inline std::size_t hash_value(const atom_type2_tuple& tup)
{

	/* old version
	boost::hash<char[(atom_type::MAX_LABEL_LEN * 2)]> hasher;
	char str[(atom_type::MAX_LABEL_LEN * 2)];
	for (unsigned int i = 0; i < atom_type::MAX_LABEL_LEN; i++){
		str[i] = tup.get<0>().atom_label[i];
	}
	for (unsigned int i = 0; i < atom_type::MAX_LABEL_LEN; i++){
		str[i+atom_type::MAX_LABEL_LEN] = tup.get<1>().atom_label[i];
	}

	return hasher(str);
	*/

	std::size_t seed = 0;

	boost::hash_combine(seed, tup.get<0>());
	boost::hash_combine(seed, tup.get<1>());

	return seed;

}

inline atom_type::atom_type(){
	for (unsigned int i = 0; i < MAX_LABEL_LEN; i++){
		this->atom_label[i] = ' ';
	}
	this->set_ff_atom_id(UNK_TYPE);
}

inline atom_type::atom_type(const std::string inlabel, const bool noff){
	std::string label = inlabel;
	boost::trim(label);
	const unsigned  int str_len = label.length();
	for (unsigned int i = 0; i < MAX_LABEL_LEN; i++){
		if (i < str_len){
			this->atom_label[i] = label[i];
		}
		else {
			this->atom_label[i] = ' ';
		}
	}
	if (!noff){
		this->auto_set_ff_atom_id();
	}
	else {
		this->set_ff_atom_id(UNK_TYPE);
	}
}

inline atom_type& atom_type::operator=(const atom_type& at){
	for (unsigned int i = 0; i < MAX_LABEL_LEN; i++){
		this->atom_label[i] = at.atom_label[i];
	}
	this->ff_atom_id = at.ff_atom_id;
	return *this;
}

inline std::string atom_type::get_label() const{
	std::string ret_str(MAX_LABEL_LEN, ' ');
	for (unsigned int i = 0; i < MAX_LABEL_LEN; i++){
		ret_str[i] = this->atom_label[i];
	}
	return ret_str;
}

inline std::string atom_type::get_pdb_formatted_label() const{
	std::string ret_str = this->get_label();
	boost::trim(ret_str);
	unsigned int size = ret_str.size();

	if (size == MAX_LABEL_LEN){
		return ret_str;
	}
	else if (size == 0){
		return std::string("    ");
	}
	else {
		if (std::isdigit(ret_str[0])){
			for (unsigned int i = size; i < MAX_LABEL_LEN; i++){
				ret_str.append(" ");
			}
			return ret_str;
		}
		else {
			ret_str.insert(0," ");
			size = ret_str.size();
			for (unsigned int i = size; i < MAX_LABEL_LEN; i++){
				ret_str.append(" ");
			}
			return ret_str;
		}
	}
	return ret_str;
}

inline void atom_type::set_ff_atom_id(const int val){
	this->ff_atom_id = val;
}

inline int atom_type::get_ff_atom_id() const{
	return this->ff_atom_id;
}

//! holding function for now
/*
inline std::string atom_type::get_ff_string() const{
	return this->get_label();
}
*/

inline bool operator==( const atom_type& at1,  const atom_type& at2 ){
	for (unsigned int i = 0; i < atom_type::MAX_LABEL_LEN; i++){
		if (at1.atom_label[i] != at2.atom_label[i]){
			return false;
		}
	}
	return true;
}

inline bool operator!=( const atom_type& at1,  const atom_type& at2){
	return !(at1 == at2);
}

inline bool operator<( const atom_type& at1,  const atom_type& at2 ){
	for (unsigned int i = 0; i < atom_type::MAX_LABEL_LEN; i++){
		if (at1.atom_label[i] < at2.atom_label[i]){
			return true;
		}
		else if (at1.atom_label[i] > at2.atom_label[i]) {
			return false;
		}
	}
	return false;
}



}
}


#endif /* ATOM_TYPE_H_ */
