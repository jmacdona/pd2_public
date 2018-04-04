/*
 * ca_pseudo_cb_pos.h
 *
 *  Created on: 22 Feb 2011
 *      Author: jmacdona
 */

#ifndef CA_PSEUDO_CB_POS_H_
#define CA_PSEUDO_CB_POS_H_


#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"
#include <map>


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

class ca_pseudo_cb_pos;

typedef boost::shared_ptr<ca_pseudo_cb_pos> ca_pseudo_cb_pos_shared_ptr;

ca_pseudo_cb_pos_shared_ptr new_ca_pseudo_cb_pos();

class ca_pseudo_cb_pos : public potential_interface {


	friend ca_pseudo_cb_pos_shared_ptr new_ca_pseudo_cb_pos();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool init();

	void add_rst(const int resnum, const UTILS::vector3d vec);

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

private:
	ca_pseudo_cb_pos();
	std::map<int, UTILS::vector3d> resnum_ca_pos_map;




};


}
}
}
}




#endif /* CA_PSEUDO_CB_POS_H_ */
