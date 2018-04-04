/*
 * ca_GO.h
 *
 *  Created on: 3 Oct 2013
 *      Author: jmacdona
 */

#ifndef CA_GO_H_
#define CA_GO_H_
#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{


potential_shared_ptr new_ca_GO();

class ca_GO : public potential_interface {



	friend potential_shared_ptr new_ca_GO();

public:

	virtual ~ca_GO(){
	}

	ca_GO();
	double upper_limit;
	int min_seq_sep;

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

};

}
}
}
}



#endif /* CA_GO_H_ */
