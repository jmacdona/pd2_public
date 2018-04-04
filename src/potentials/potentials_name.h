//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * potentials_name.h
 *
 *  Created on: 17 Mar 2010
 *      Author: jmacdona
 */

#ifndef POTENTIALS_NAME_H_
#define POTENTIALS_NAME_H_
#include <string>
#include <cctype>
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
#include <map>
#include <vector>

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

const int max_label_size = 20;

//! stores potential name which must be less than or equal to 20 characters long
class potentials_name {




	//! comparison using potentials_label[]
	friend bool operator==( const potentials_name&,  const potentials_name& );
	//! comparison using potentials_label[]
	friend bool operator<( const potentials_name&,  const potentials_name& );

public:

	//potentials_name();
	potentials_name(const std::string label);
	potentials_name& operator=(const potentials_name&);


	std::string get_label() const;


private:

	void set_label(const std::string& label);

	//! name as read by weights so must be unique and less than 20 chars long.
	char potentials_label[max_label_size+1];

};

typedef std::map<potentials_name, double> potentials_name_double_map;
typedef std::vector<potentials_name> potentials_name_vector;


inline potentials_name::potentials_name(const std::string label){

	this->set_label(label);

}

inline potentials_name& potentials_name::operator=(const potentials_name& at){
	for (int i = 0; i < max_label_size; i++){
		this->potentials_label[i] = at.potentials_label[i];
	}

	return *this;
}

inline bool operator==( const potentials_name& at1,  const potentials_name& at2 ){
	for (int i = 0; i < max_label_size; i++){
		if (at1.potentials_label[i] != at2.potentials_label[i]){
			return false;
		}
	}
	return true;
}


inline bool operator<( const potentials_name& at1,  const potentials_name& at2 ){
	for (int i = 0; i < max_label_size; i++){
		if (at1.potentials_label[i] < at2.potentials_label[i]){
			return true;
		}
		else if (at1.potentials_label[i] > at2.potentials_label[i]) {
			return false;
		}
	}
	return false;
}


}
}
}

#endif /* POTENTIALS_NAME_H_ */
