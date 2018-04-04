/*
 * simple_ca_mc_protocol.h
 *
 *  Created on: 23 Aug 2010
 *      Author: jmacdona
 */

#ifndef SIMPLE_CA_MC_PROTOCOL_H_
#define SIMPLE_CA_MC_PROTOCOL_H_
#include "simulation/mc_protocol_interface.h"
#include <limits>

namespace PRODART {
namespace POSE {
namespace SIM {
namespace CA {


class simple_ca_mc_protocol;


typedef boost::shared_ptr<simple_ca_mc_protocol> simple_ca_mc_protocol_shared_ptr;

mc_protocol_interface_shared_ptr new_simple_ca_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr);

class simple_ca_mc_protocol : public mc_protocol_interface {


	friend mc_protocol_interface_shared_ptr new_simple_ca_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr);

private:


protected:


	simple_ca_mc_protocol(unsigned long steps, double _beta, MTRand::MTRand_shared_ptr mt_ptr) : mc_protocol_interface(steps, _beta, mt_ptr){
		/*
		this->steps_to_run = steps;
		this->beta = _beta;
		*/
	}

public:




};





}
}
}
}

#endif /* SIMPLE_CA_MC_PROTOCOL_H_ */
