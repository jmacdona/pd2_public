/*
 * lj_nb_potential.cpp
 *
 *  Created on: 22 Nov 2010
 *      Author: jmacdona
 */

#include "lj_nb_potential.h"

using namespace boost;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;
using namespace std;

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

const double lj_nb_potential::scale_to_joules = 4184; //scale kcal to J
const double lj_nb_potential::R = 8.314472; //J K-1 mol-1
const double lj_nb_potential::assumed_temperature = 200.0;
const double lj_nb_potential::scale_factor = (scale_to_joules / (R * assumed_temperature));


potential_shared_ptr new_lj_nb_potential(){
	potential_shared_ptr ptr(new lj_nb_potential());
	return ptr;
}

lj_nb_potential::lj_nb_potential() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("lj_nb_potential"));
}




bool lj_nb_potential::init(){

	return true;
}


}
}
}


