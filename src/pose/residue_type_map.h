//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * residue_type_map.h
 *
 *  Created on: Feb 20, 2010
 *      Author: jmacdon
 */

#ifndef RESIDUE_TYPE_MAP_H_
#define RESIDUE_TYPE_MAP_H_


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>
#include <vector>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/thread.hpp>

#include "residue_type.h"

#include "prodart_env/prodart_env.h"


namespace PRODART {
namespace POSE {

enum StandardResidueType {UNDEF_RES, UNK_RES, ALA, VAL, PHE, PRO, MET, ILE, LEU, ASP, GLU, LYS, ARG, SER,
	THR, TYR, HIS, CYS, ASN, GLN, TRP, GLY};

typedef std::map<residue_type, char> residue_type_char_map;
typedef std::map<char, residue_type> char_residue_type_map;


//! singleton class
class residue_type_map{

private:
	residue_type_map(const residue_type_map&);
	residue_type_char_map one_letter_name_map;
	char_residue_type_map char_to_res_type_map;

	static void Init();



protected:
	residue_type_map();
	virtual ~residue_type_map();

public:

	static const residue_type_map* Instance();

	char to_1_letter_code(const residue_type& rt) const;
	std::string to_1_letter_code(const residue_type_vector& rt_vec) const;
	residue_type to_residue_type(char one_letter_code) const;

	bool is_standard_type(const residue_type& ty) const;


};

inline char residue_type_map::to_1_letter_code(const residue_type& rt) const{
	//return one_letter_name_map[rt];
	return one_letter_name_map.find(rt) != one_letter_name_map.end() ? one_letter_name_map.find(rt)->second : 'U' ;

}

inline std::string residue_type_map::to_1_letter_code(const residue_type_vector& rt_vec) const{
	//return one_letter_name_map[rt];
	std::string ret_str;
	residue_type_vector::const_iterator it;
	for (it = rt_vec.begin(); it != rt_vec.end(); it++){
		ret_str.push_back(to_1_letter_code(*it));
	}
	return ret_str;
}

inline residue_type residue_type_map::to_residue_type(char one_letter_code) const{
	return char_to_res_type_map.find(one_letter_code) != char_to_res_type_map.end() ? char_to_res_type_map.find(one_letter_code)->second : residue_type("UNK") ;
}

inline bool residue_type_map::is_standard_type(const residue_type& ty) const{
	if (ty == residue_type("UNK") || ty == residue_type("") || ty == residue_type("ASX") ){
		return false;
	}
	else {
		if (one_letter_name_map.find(ty) != one_letter_name_map.end()){
			return true;
		}
		else {
			return false;
		}
	}
}

residue_type_vector one_letter_seq_to_residue_type_vector(std::string seq);


}
}

#endif /* RESIDUE_TYPE_MAP_H_ */
