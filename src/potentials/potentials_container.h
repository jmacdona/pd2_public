//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * potentials_container.h
 *
 *  Created on: 11 Jun 2010
 *      Author: jmacdona
 */

#ifndef POTENTIALS_CONTAINER_H_
#define POTENTIALS_CONTAINER_H_


#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

class potentials_container;

typedef boost::shared_ptr<potentials_container> potentials_container_shared_ptr;

potentials_container_shared_ptr new_potentials_container();

class potentials_container : public potential_interface {

	friend potentials_container_shared_ptr new_potentials_container();
	friend class potentials_factory;

	potentials_energies_map default_energies_map;

protected:

	potentials_container();

	potential_shared_ptr_vector potentials;


public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	void get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			bool_vector& vec) const;

	double get_energy_residue_loop(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
				potentials_energies_map& energies_map,
				const bool_vector& res_loop_mask) const;

	bool init();


	//! does it provide energies labeled with this name
	bool provides(const potentials_name& query_name) const;

	//! add potential to container but will not have corresponding component weights entries in the energies_map so you'll need to add this manually
	void add_potential(potential_shared_ptr ptr);

	void add_potential(potential_shared_ptr ptr, const double weight);

	potential_shared_ptr get_potential(const potentials_name& name);

	const potentials_energies_map& get_default_energies_map() const{
		return default_energies_map;
	}

	potentials_energies_map& get_default_energies_map(){
		return default_energies_map;
	}

private:



};




}
}
}

#endif /* POTENTIALS_CONTAINER_H_ */
