/*
 * ca_frag_freq.h
 *
 *  Created on: 27 Aug 2010
 *      Author: jmacdona
 */

#ifndef CA_FRAG_FREQ_H_
#define CA_FRAG_FREQ_H_

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

potential_shared_ptr new_ca_frag_freq();

class ca_frag_freq : public potential_interface {


	friend potential_shared_ptr new_ca_frag_freq();

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
	ca_frag_freq();
	//const PRODART::POSE::POTENTIALS::potentials_name name;

	double get_frag_energy(const int prev_prev_frag_num, const int prev_frag_num, const int frag_num) const;

	int get_combined_bin(const int n1, const int n2) const;
	void deconvolute_combined_bin(const int bin,int& n1, int& n2) const;

	int get_combined_bin(const int n1, const int n2, const int n3) const;
	void deconvolute_combined_bin(const int bin,int& n1, int& n2, int& n3) const;

	//int_double_map frag_energy;
	//int_double_map frag12_energy;
	//int_double_map frag13_energy;
	//int_double_map frag123_energy;

	int_double_umap frag_energy;
	int_double_umap frag12_energy;
	int_double_umap frag13_energy;
	int_double_umap frag123_energy;

	int total_num_frags;



};


}
}
}
}










#endif /* CA_FRAG_FREQ_H_ */
