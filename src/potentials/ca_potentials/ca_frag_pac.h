/*
 * ca_frag_pac.h
 *
 *  Created on: Feb 27, 2011
 *      Author: jmacdon
 */

#ifndef CA_FRAG_PAC_H_
#define CA_FRAG_PAC_H_




#include "../potential_interface.h"
#include "pose_meta/ca_pose_meta.h"


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

//typedef std::vector< double > double_vector;
//typedef std::map< int, double_vector, std::less< int > > int_double_vector_map;
typedef std::map< int, double, std::less< int > > int_double_map;
typedef std::map< POSE::residue_type, int > residue_type_int_map;

potential_shared_ptr new_ca_frag_pac();

class ca_frag_pac : public potential_interface {


	friend potential_shared_ptr new_ca_frag_pac();

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
			potentials_energies_map& energies_map) const;

	bool init();

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;

	std::istream& load_data( std::istream& input );

protected:

	residue_type_int_map restype_to_int_map;
	int_double_map PAC_score_map;



	ca_frag_pac();

	int rt_to_int(const POSE::residue_type &rt) const;
	int overall_bin(const int frag_num,
				const POSE::residue_type &rt1,
				const POSE::residue_type &rt2) const;

};


}
}
}
}








#endif /* CA_FRAG_PAC_H_ */
