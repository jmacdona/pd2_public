/*
 * minimiser_interface.h
 *
 *  Created on: 15 Sep 2013
 *      Author: jmacdona
 */

#ifndef MINIMISER_INTERFACE_H_
#define MINIMISER_INTERFACE_H_

#include "prodart_env/prodart_env.h"
#include "movers/mover_interface.h"
#include "pose_meta/pose_meta_interface.h"
#include <boost/shared_ptr.hpp>

namespace PRODART {
namespace POSE {
namespace SIM {


class minimiser_interface;
typedef boost::shared_ptr<minimiser_interface> minimiser_interface_shared_ptr;

class minimiser_interface : public PRODART::POSE::MOVERS::mover_interface {

private:



protected:

	minimiser_interface(){
	}

	virtual ~minimiser_interface(){
	}

public:

	virtual bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const = 0;

	virtual void set_steps(unsigned int val) = 0;
	virtual void set_tol(double tol) = 0;
	virtual void set_lin_min_tol(double tol) = 0;

};




}
}
}



#endif /* MINIMISER_INTERFACE_H_ */
