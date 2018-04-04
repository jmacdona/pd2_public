/*
 * simple_ca_sim_anneal_protocol.cpp
 *
 *  Created on: 15 Oct 2010
 *      Author: jmacdona
 */

#include "simple_ca_sim_anneal_protocol.h"



namespace PRODART {
namespace POSE {
namespace SIM {
namespace CA {



mc_protocol_interface_shared_ptr new_simple_ca_sim_anneal_protocol(unsigned long start_beta_steps_,
		unsigned long anneal_steps_,
		unsigned long final_beta_steps_,
		double start_beta_,
		double final_beta_,
		MTRand::MTRand_shared_ptr mt_ptr){
	mc_protocol_interface_shared_ptr ptr(new simple_ca_sim_anneal_protocol(start_beta_steps_,
			anneal_steps_,
			final_beta_steps_,
			start_beta_,
			final_beta_,
			mt_ptr));
	return ptr;
}






}
}
}
}

