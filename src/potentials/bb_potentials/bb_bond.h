/*
 * bb_bond.h
 *
 *  Created on: 25 Oct 2010
 *      Author: jmacdona
 */

#ifndef BB_BOND_H_
#define BB_BOND_H_

#include "potentials/bond_potential.h"
#include "pose_meta/bb_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

potential_shared_ptr new_bb_bond();

class bb_bond : public POTENTIALS::bond_potential {

	friend potential_shared_ptr new_bb_bond();

protected:

	bb_bond();


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



#endif /* BB_BOND_H_ */
