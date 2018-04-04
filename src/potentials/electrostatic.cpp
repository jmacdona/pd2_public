/*
 * electrostatic.cpp
 *
 *  Created on: 23 Nov 2010
 *      Author: jmacdona
 */
#include "electrostatic.h"

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

const double electrostatic::scale_to_joules = 4184; //scale kcal to J
const double electrostatic::R = 8.314472; //J K-1 mol-1
const double electrostatic::assumed_temperature = 200.0;
const double electrostatic::scale_factor = (electrostatic::scale_to_joules / (electrostatic::R * electrostatic::assumed_temperature));

const double electrostatic::e_r = 5.0;		//80.0 in water
const double electrostatic::e0 = 8.854187817e-12; //C2 N−1 m−2 or C2/(J.m)
const double electrostatic::avogadro = 6.0221415e23;
const double electrostatic::elementary_charge = 1.602176487e-19;
const double electrostatic::k_ele = (electrostatic::elementary_charge*electrostatic::elementary_charge*electrostatic::avogadro * (1e10)
	* (1.0/electrostatic::scale_to_joules)
	*(1.0/(4.0*PI*electrostatic::e0))) / electrostatic::e_r;		//  = 332 to get kcal mol-1
const double electrostatic::overall_ele_constant = electrostatic::k_ele / electrostatic::e_r;


potential_shared_ptr new_electrostatic(){
	potential_shared_ptr ptr(new electrostatic());
	return ptr;
}

electrostatic::electrostatic() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("electrostatic"));
}




bool electrostatic::init(){

	return true;
}


}
}
}


