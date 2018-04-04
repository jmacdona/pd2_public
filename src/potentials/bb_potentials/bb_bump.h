/*
 * bb_bump.h
 *
 *  Created on: 26 Oct 2010
 *      Author: jmacdona
 */

#ifndef BB_BUMP_H_
#define BB_BUMP_H_

#include "pose_meta/bb_pose_meta.h"
#include "../potential_interface.h"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <map>
#include <boost/unordered_map.hpp>
#include <cstring>

namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace BB{

typedef std::map<POSE::atom_type, double> atom_type_double_map;
typedef boost::tuple<POSE::atom_type, POSE::atom_type> atom_type2_tuple;
typedef std::map<atom_type2_tuple, double> atom_type2_tuple_double_map;

/*
struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};
*/

typedef boost::unordered_map<std::string, double> string_double_unord_map;
typedef boost::unordered_map<atom_type2_tuple, double, boost::hash<PRODART::POSE::atom_type2_tuple> > atom_type2_tuple_double_unord_map;
typedef boost::unordered_map<POSE::atom_type, double, boost::hash<PRODART::POSE::atom_type> > atom_type_double_unord_map;

potential_shared_ptr new_bb_bump();

class bb_bump : public POTENTIALS::potential_interface {

	friend potential_shared_ptr new_bb_bump();

protected:

	bb_bump();

	double get_energy(const PRODART::POSE::META::nb_pair_element& ele) const;
	double get_energy_grad(const PRODART::POSE::META::nb_pair_element& ele, UTILS::vector3d_vector& grad, const double weight ) const;

	double get_min_dist(const POSE::atom_type, const POSE::atom_type) const;
	double get_min_dist_sq(const POSE::atom_type, const POSE::atom_type) const;

	void add_min_dist(const POSE::atom_type at1, const POSE::atom_type at2, const double dist );


	const int min_bond_sep;
	const double vdw_deflator;
	const double default_min_dist, default_min_dist_sq;
	const double outer_bound;

	//atom_type_double_map atomVdwParams;
	atom_type_double_unord_map atomVdwParams;

	//atom_type2_tuple_double_map min_dists_map;
	//atom_type2_tuple_double_map min_dists_sq_map;

	atom_type2_tuple_double_unord_map min_dists_map;
	atom_type2_tuple_double_unord_map min_dists_sq_map;

	string_double_unord_map min_dists_un_map;
	string_double_unord_map min_dists_sq_un_map;

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


#endif /* BB_BUMP_H_ */
