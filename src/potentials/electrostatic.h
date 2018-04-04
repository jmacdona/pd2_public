/*
 * electrostatic.h
 *
 *  Created on: 23 Nov 2010
 *      Author: jmacdona
 */

#ifndef ELECTROSTATIC_H_
#define ELECTROSTATIC_H_


#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <map>

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

typedef boost::tuple<double, double> double2_tuple;
typedef std::map<POSE::atom_type, double2_tuple> atom_type_double2_tuple_map;
typedef boost::tuple<POSE::atom_type, POSE::atom_type> atom_type2_tuple;
typedef std::map<atom_type2_tuple, double2_tuple> atom_type2_tuple_double2_tuple_map;

potential_shared_ptr new_electrostatic();

class electrostatic : public potential_interface {


	friend potential_shared_ptr new_electrostatic();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const{
		return 0;
	}

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const{
		return 0;
	}


	bool init();

	static const double scale_to_joules;// = 4184; //scale kcal to J
	static const double R;// = 8.314472; //J K-1 mol-1
	static const double assumed_temperature;// = 200.0;
	static const double scale_factor;// = (scale_to_joules / (R * assumed_temperature));

	static const double e_r ;//= 5.0;		//80.0 in water
	static const double e0;//=8.854187817e-12; //C2 N−1 m−2 or C2/(J.m)
	static const double avogadro;// = 6.0221415e23;
	static const double elementary_charge;// = 1.602176487e-19;
	static const double k_ele;// = (elementary_charge*elementary_charge*avogadro * (1e10)
		//* (1.0/scale_to_joules)
		//*(1.0/(4.0*PI*e0))) / e_r;		//  = 332 to get kcal mol-1
	static const double overall_ele_constant;// = geometricPot::k_ele / geometricPot::e_r;


protected:
	electrostatic();




};


}
}
}





#endif /* ELECTROSTATIC_H_ */
