/*
 * fixed_bb_motif_anneal.cpp
 *
 *  Created on: Mar 25, 2012
 *      Author: jmacdona
 */

#include "fixed_bb_motif_anneal.h"

namespace PRODART {
namespace POSE {
namespace SIM {
namespace BB {



mc_protocol_interface_shared_ptr new_fixed_bb_motif_anneal(unsigned long start_beta_steps_,
		unsigned long anneal_steps_,
		unsigned long final_beta_steps_,
		double start_beta_,
		double final_beta_,
		MTRand::MTRand_shared_ptr mt_ptr){
	mc_protocol_interface_shared_ptr ptr(new fixed_bb_motif_anneal(start_beta_steps_,
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


