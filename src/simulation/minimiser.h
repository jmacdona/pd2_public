/*
 * minimiser.h
 *
 *  Created on: 28 Sep 2010
 *      Author: jmacdona
 */

#ifndef MINIMISER_H_
#define MINIMISER_H_

#ifdef PD2_USE_CERES
#include "ceresminimiser.h"
#else
#include "gslminimiser.h"
#endif


namespace PRODART {
	namespace POSE {
		namespace SIM {
			
			
			//return default minimiser
			minimiser_interface_shared_ptr new_minimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use);
			//return default minimiser
			minimiser_interface_shared_ptr new_minimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection);
			
			
		}
	}
}



#endif /* MINIMISER_H_ */
