/*
 * bb_14lj.h
 *
 *  Created on: 22 Nov 2010
 *      Author: jmacdona
 */

#ifndef BB_14LJ_H_
#define BB_14LJ_H_

#include "potentials/lj_nb_potential.h"
#include "pose_meta/bb_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

potential_shared_ptr new_bb_14_15lj();

class bb_14_15lj : public POTENTIALS::lj_nb_potential {

	friend potential_shared_ptr new_bb_14_15lj();

protected:

	bb_14_15lj();


public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

private:



};

}
}
}
}



#endif /* BB_14LJ_H_ */
