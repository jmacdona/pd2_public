//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * potentials_store.h
 *
 *  Created on: 17 Mar 2010
 *      Author: jmacdona
 */

#ifndef POTENTIALS_STORE_H_
#define POTENTIALS_STORE_H_

#include <string>
#include <cctype>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>
#include "potentials_name.h"
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>




namespace PRODART {
namespace POSE {
namespace POTENTIALS {

class potentials_store {

public:

	double weight;
	double energy;


	potentials_store() : weight(0), energy(0){

	}
	potentials_store(double weight_, double energy_) : weight(weight_), energy(energy_){

	}

	potentials_store operator-(const potentials_store& st2) const {
		potentials_store rtn_st(this->weight - st2.weight, this->energy - st2.energy);
		return rtn_st;
	}

};

typedef std::map<potentials_name, potentials_store> potentials_name_potentials_store_map;//potentials_energies_map;

class potentials_energies_map : public  potentials_name_potentials_store_map{

public:

	double total_energy;

	potentials_energies_map() : total_energy(0){

	}

	bool double_isOK( const double val, const potentials_name label ) const{
		bool isOK = true;
		if (boost::math::isnan(val)){
			std::cerr << "add_energy_component: ERROR: isnan at label: " << label.get_label() << std::endl;
			isOK = false;
		}
		else if (boost::math::isinf(val)){
			std::cerr << "add_energy_component: ERROR: isinf at label: " << label.get_label() << std::endl;
			isOK = false;
		}
		return isOK;
	}

	double add_energy_component(potentials_name label, const double add_energy){

		potentials_store& store = (*this)[label];
		assert(double_isOK(add_energy, label));

		store.energy = add_energy;
		const double wt_energy = store.weight * store.energy;
		this->total_energy += wt_energy;

		assert(double_isOK(wt_energy, label));

		return wt_energy;
	}

	double get_weight(potentials_name label) const{

		if (this->find(label) != this->end()){
			const potentials_store& store = this->find(label)->second;
			return store.weight;
		}
		return 0;
	}

	void adjust_weight(const potentials_name label, const double weight) {
		if (this->find(label) != this->end()){
			potentials_store& store = this->find(label)->second;
			store.weight = weight;
		}
		else {
			std::cerr << "potentials_energies_map: ERROR: weight not in map:\t" << label.get_label() << std::endl;
		}
		return;
	}

	bool contains_term(potentials_name label) const{
		if (this->find(label) != this->end()){
			return true;
		}
		return false;
	}

	potentials_energies_map get_diff(const potentials_energies_map other) const;


	std::ostream& print_headers(std::ostream& output) const;
	std::ostream& print_unweighted_components(std::ostream& output) const;
	std::ostream& print_weighted_components(std::ostream& output) const;
	std::ostream& print_weights(std::ostream& output) const;

	std::ostream& print_unweighted_diff(std::ostream& output, const potentials_energies_map other) const;


};


}
}
}





#endif /* POTENTIALS_STORE_H_ */
