/*
 * simple_ca_mc_protocol.cpp
 *
 *  Created on: 23 Aug 2010
 *      Author: jmacdona
 */

#include "simple_ca_mc_protocol.h"


namespace PRODART {
namespace POSE {
namespace SIM {
namespace CA {



mc_protocol_interface_shared_ptr new_simple_ca_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr){
	mc_protocol_interface_shared_ptr ptr(new simple_ca_mc_protocol(steps, _beta, mt_ptr));
	return ptr;
}






}
}
}
}

