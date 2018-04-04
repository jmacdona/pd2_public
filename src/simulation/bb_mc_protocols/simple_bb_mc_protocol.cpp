/*
 * simple_bb_mc_protocol.cpp
 *
 *  Created on: 8 Nov 2010
 *      Author: jmacdona
 */
#include "simple_bb_mc_protocol.h"


namespace PRODART {
namespace POSE {
namespace SIM {
namespace BB {



mc_protocol_interface_shared_ptr new_simple_bb_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr){
	mc_protocol_interface_shared_ptr ptr(new simple_bb_mc_protocol(steps, _beta, mt_ptr));
	return ptr;
}






}
}
}
}


