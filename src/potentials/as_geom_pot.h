/*
 * as_geom_pot.h
 *
 *  Created on: Jan 2, 2011
 *      Author: jmacdon
 */

#ifndef AS_GEOM_POT_H_
#define AS_GEOM_POT_H_

#include "potential_interface.h"
#include "pose_meta/pose_meta_interface.h"
#include "pose_meta/ca_pose_meta.h"
#include "pose_meta/bb_pose_meta.h"
#include "pose_utils/pose_utils.h"
#include <vector>
#include <map>
#include <boost/tuple/tuple.hpp>
#include <sstream>
#include <numeric>
#include "utils/math_utils.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

typedef std::vector< double > double_vector;
typedef std::vector< int > int_vector;
typedef std::map<int, int> int_int_map;
typedef std::vector< int_int_map > int_int_map_vector;
typedef boost::tuple<std::string, char, std::string, char> string_char_string_char_tuple;
typedef std::vector<string_char_string_char_tuple> string_char_string_char_tuple_vector;

class as_geom_pot;

typedef boost::shared_ptr<as_geom_pot> as_geom_pot_shared_ptr;

potential_shared_ptr new_as_geom_pot();

class as_geom_pot : public potential_interface {


	friend potential_shared_ptr new_as_geom_pot();

protected:

	as_geom_pot();

	void temp_store_mapping();
	void temp_restore_mapping();


	void dummy_funct(double val) const;

	bool motif_loaded;

	POSE::pose_shared_ptr motif_pdb;

	static const int NO_MAPPING_ = 9999999;
	static const double grad_h;

	int_int_map_vector motif_target_res_map_vector;
	int_int_map_vector temp_bk_motif_target_res_map_vector;
	int_int_map_vector move_bk_motif_target_res_map_vector;

	string_char_string_char_tuple_vector fixed_mapping;

	MTRand::MTRand_shared_ptr rand_num;

	bool_vector loop_mask;

	bool using_fixed_mapping;

	const prot_backbone_map* const bb_map;// bb_map(prot_backbone_map::Instance());

	bool assign_motif_by_first_res(POSE::pose_shared_ptr pose_, const int_vector& vec);

public:

	double get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;

	double get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
			potentials_energies_map& energies_map) const;


	bool init();

	//! load motif pdb
	std::istream& load_data( std::istream& input );
	std::istream& load_mapping( std::istream& input );

	bool make_move(POSE::pose_shared_ptr pose_);

	bool move_element(POSE::pose_shared_ptr pose_,
			const int ele_num,
			const int max_move);

	bool apply_fixed_mapping(POSE::pose_shared_ptr pose_);


	void move_store_mapping();
	void move_restore_mapping();

	bool random_assign_motif(POSE::pose_shared_ptr pose_);
	void print_assign_info(std::ostream& output) const;

	void get_ca_atom_mapping(POSE::pose_shared_ptr pose_,
			std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> &atom_mapping,
			pose_shared_ptr this_motif = pose_shared_ptr()) const;
	double get_ca_rmsd(POSE::pose_shared_ptr pose_) const;

	void get_bb_atom_mapping(POSE::pose_shared_ptr pose_,
			std::map<PRODART::POSE::atom_shared_ptr, PRODART::POSE::atom_shared_ptr> &atom_mapping,
			pose_shared_ptr this_motif = pose_shared_ptr()) const;
	double get_bb_rmsd(POSE::pose_shared_ptr pose_) const;

	bool is_loaded() const{
		return  motif_loaded;
	}

	const int_int_map_vector& get_motif_assignment() const{
		return motif_target_res_map_vector;
	}

	void set_motif_assignment(const int_int_map_vector& val){
		motif_target_res_map_vector = val;
	}

	void set_loop_mask(const bool_vector& vec){
		loop_mask.clear();
		loop_mask.resize(vec.size(), false);
		for (unsigned int i = 0; i < vec.size(); i++){
			loop_mask[i] = !vec[i];
		}
	}

	pose_shared_ptr get_aligned_motif(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_) const;

	void paint_motif_occupancy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_, const double base_value = 1.0 ) ;

	//! does it provide energies labeled with this name
	//bool provides(const potentials_name& query_name) const;


	//! copy sidechains from motif to pose_
	void transfer_sidechains(const PRODART::POSE::META::bb_pose_meta_shared_ptr pose_meta_) const;
	void transfer_non_peptides(const PRODART::POSE::META::bb_pose_meta_shared_ptr pose_meta_) const;

	//returns RMSD
	double do_exaustive_bb_motif_search(POSE::pose_shared_ptr pose_);


};

}
}
}


#endif /* AS_GEOM_POT_H_ */
