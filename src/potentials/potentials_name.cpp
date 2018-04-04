//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * potentials_name.cpp
 *
 *  Created on: 17 Mar 2010
 *      Author: jmacdona
 */

#include "potentials_name.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {






void potentials_name::set_label(const std::string& label){

	const int str_len = label.length();
	for (int i = 0; i < max_label_size; i++){
		if (i < str_len){
			this->potentials_label[i] = label[i];
		}
		else {
			this->potentials_label[i] = '\0';
		}
	}
	this->potentials_label[max_label_size] = '\0';

}


std::string potentials_name::get_label() const{
	return std::string(potentials_label);
}


}
}
}

