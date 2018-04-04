/*
 * fragment_classifier_interface.h
 *
 *  Created on: 25 Aug 2010
 *      Author: jmacdona
 */

#ifndef FRAGMENT_CLASSIFIER_INTERFACE_H_
#define FRAGMENT_CLASSIFIER_INTERFACE_H_
#include "pose_meta/pose_meta_interface.h"
#include "pose_meta/pose_meta_defs.h"
#include "pose_meta/ca_pose_meta.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include "prodart_env/prodart_env.h"
#include <exception>

namespace PRODART {
namespace POSE {
namespace META {
namespace FRAG_CLASS {

class fragment_classifier_interface;

typedef boost::shared_ptr<fragment_classifier_interface> fragment_classifier_shared_ptr;

class fragment_classifier_interface {

public:
	virtual bool classify_fragments (const pose_shared_ptr _pose,
			const PRODART::POSE::META::ca_pose_meta_shared_ptr _pose_meta) const = 0;

protected:
	fragment_classifier_interface(){

	}
	virtual ~fragment_classifier_interface(){
	}

};



}
}
}
}

#endif /* FRAGMENT_CLASSIFIER_INTERFACE_H_ */
