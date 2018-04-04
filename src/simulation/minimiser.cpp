/*
 * minimiser.cpp
 *
 *  Created on: 28 Sep 2010
 *      Author: jmacdona
 */

#include "minimiser.h"



using namespace PRODART;
using namespace PRODART::UTILS;
using namespace PRODART::POSE;

using namespace std;


namespace PRODART {
namespace POSE {
namespace SIM {

	//return default minimiser
	minimiser_interface_shared_ptr new_minimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use){
		//return new_gslminimiser(potentials_to_use);
#ifdef PD2_USE_CERES
		return new_ceresminimiser(potentials_to_use);
#else
		return new_gslminimiser(potentials_to_use);
#endif
	}
	//return default minimiser
	minimiser_interface_shared_ptr new_minimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection){
		//return new_gslminimiser(potentials_to_use, atom_selection);
#ifdef PD2_USE_CERES
		return new_ceresminimiser(potentials_to_use, atom_selection);
#else
	return new_gslminimiser(potentials_to_use, atom_selection);	
#endif

	}

}
}
}



