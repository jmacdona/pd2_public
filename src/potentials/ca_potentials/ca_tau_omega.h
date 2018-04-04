/*
 * ca_tau_omega.h
 *
 *  Created on: 26 Aug 2010
 *      Author: jmacdona
 */

#ifndef CA_TAU_OMEGA_H_
#define CA_TAU_OMEGA_H_

#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"
#include <boost/unordered_map.hpp>


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

typedef std::vector< double > double_vector;
typedef std::map< int, double_vector, std::less< int > > int_double_vector_map;
typedef std::map< int, double, std::less< int > > int_double_map;

typedef boost::unordered_map< int, double > int_double_umap;

potential_shared_ptr new_ca_tau_omega();

class ca_tau_omega : public potential_interface {


	friend potential_shared_ptr new_ca_tau_omega();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

	std::istream& load_data( std::istream& input );

private:
	ca_tau_omega();
	//const PRODART::POSE::POTENTIALS::potentials_name name;

	double frag_energy(const PRODART::POSE::META::frag4_element& ele) const;
	double get_dihedral_distance(const double dih1, const double dih2) const;

	int_double_umap tau_0_params;
	int_double_umap tau_k_params;
	int_double_umap tau_const_params;

	int_double_umap omega1_0_params;
	int_double_umap omega1_k_params;
	int_double_umap omega1_const_params;

	int_double_umap omega2_0_params;
	int_double_umap omega2_k_params;
	int_double_umap omega2_const_params;


};


}
}
}
}




#endif /* CA_TAU_OMEGA_H_ */
