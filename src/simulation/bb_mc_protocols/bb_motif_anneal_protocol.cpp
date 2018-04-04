/*
 * bb_motif_anneal_protocol.cpp
 *
 *  Created on: 5 Jan 2011
 *      Author: jmacdona
 */
#include "bb_motif_anneal_protocol.h"


namespace PRODART {
namespace POSE {
namespace SIM {
namespace BB {



mc_protocol_interface_shared_ptr new_bb_motif_anneal_protocol(unsigned long start_beta_steps_,
		unsigned long anneal_steps_,
		unsigned long final_beta_steps_,
		double start_beta_,
		double final_beta_,
		MTRand::MTRand_shared_ptr mt_ptr){
	mc_protocol_interface_shared_ptr ptr(new bb_motif_anneal_protocol(start_beta_steps_,
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


