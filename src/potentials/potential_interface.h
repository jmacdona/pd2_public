//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * potential_interface.h
 *
 *  Created on: 5 Mar 2010
 *      Author: jmacdona
 */

#ifndef POTENTIAL_INTERFACE_H_
#define POTENTIAL_INTERFACE_H_

#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "utils/angle_derivative.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "pose_meta/pose_meta_interface.h"
#include "pose_meta/ca_pose_meta.h"
#include "potentials_name.h"
#include "potentials_store.h"
#include <vector>
#include "prodart_env/prodart_env.h"
#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>


namespace PRODART {
namespace POSE {
namespace POTENTIALS {


class potential_interface;

typedef boost::shared_ptr<potential_interface> potential_shared_ptr;
typedef std::vector<potential_shared_ptr> potential_shared_ptr_vector;
typedef std::map<std::string, potential_shared_ptr> string_potential_shared_ptr_map;

typedef std::vector<std::string> string_vector;
typedef std::vector<bool> bool_vector;

class potential_interface {

protected:

	potentials_name_vector name_vector;
	bool b_is_disposable;

	potential_interface() : b_is_disposable(false){

	}

	virtual ~potential_interface(){

	}


public:


	virtual double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const = 0;

	virtual double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const = 0;

	//! for use in free energy calcs, meta-dynamics etc.
	virtual void init_set_up(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map){}

	//! for use in free energy calcs, meta-dynamics etc.
	virtual void accepted(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map){}

	virtual void get_residue_hot_spots(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			bool_vector& vec) const{
		return;
	}

	virtual double get_energy_residue_loop(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map,
			const bool_vector& res_loop_mask) const{
		return 0 ;
	}

	virtual bool training_pose(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_){
		return false;
	}

	virtual bool background_training_pose(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_){
		return false;
	}

	virtual bool calc_potentials_using_training(){
		return false;
	}

	virtual bool output_training_results(std::ostream& output){
		return false;
	}

	//! tell object to initialise - e.g. load params etc
	virtual bool init() = 0;

	//! does it provide energies labeled with this name
	virtual bool provides(const potentials_name& query_name) const {
		potentials_name_vector::const_iterator iter;
		for (iter = this->name_vector.begin(); iter != this->name_vector.end(); iter++){
			if (query_name == *iter){
				return true;
			}
		}

		return false;
	}

	//! get list of provides potential labels
	virtual potentials_name_vector provides() const{
		return name_vector;
	}

	//! tells potentials_factory to get a new instance each time
	virtual bool is_disposable() const{
		return b_is_disposable;
	}

	virtual potential_shared_ptr get_new_instance() const{
		std::cerr << "\nERROR: potential_interface: get_new_instance(): function not supported\n" << std::endl;
		return potential_shared_ptr();
	}

private:



};





}
}
}

#endif /* POTENTIAL_INTERFACE_H_ */
