/*
 * alphabet_bb_builder.h
 *
 *  Created on: May 2, 2012
 *      Author: jmacdona
 */

#ifndef ALPHABET_BB_BUILDER_H_
#define ALPHABET_BB_BUILDER_H_


#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include "prodart_env/prodart_env.h"
#include <exception>
#include "pose_meta/ca_pose_meta.h"
#include "pose_meta/frag_classify/fragment_classifier_interface.h"
#include "pose_meta/frag_classify/ca_fragment_classifier.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "pose_utils/pose_utils.h"
#include <boost/unordered_map.hpp>
#include <utility>
#include <boost/tuple/tuple.hpp>


namespace PRODART {
namespace POSE {
namespace BB_BUILDER {

class alphabet_bb_builder;

typedef boost::shared_ptr<alphabet_bb_builder> alphabet_bb_builder_shared_ptr;
typedef boost::shared_ptr<const alphabet_bb_builder> const_alphabet_bb_builder_shared_ptr;
typedef std::vector<const_pose_shared_ptr> const_pose_shared_ptr_vector;

class alphabet_bb_builder{

	friend alphabet_bb_builder_shared_ptr new_alphabet_bb_builder();


private:

	alphabet_bb_builder(const alphabet_bb_builder&);



protected:
	alphabet_bb_builder();

	const_pose_shared_ptr_vector pose_alphabet;

	int letter_length;

	double centre_weight;

	double default_b_factor;

	POSE::double_vector weights;

	bool is_even;

	double get_rmsd(const int first_res, const_pose_shared_ptr letter, pose_shared_ptr pose_) const;


	//! builds peptide bond between first_res and first_res+1. Returns rebuild atoms vector
	POSE::atom_shared_ptr_vector build_peptide_bond(const int first_res, pose_shared_ptr pose_) const;

public:
	bool load_alphabet_pose(std::istream & input);
	pose_shared_ptr get_best_fit_letter(const int first_res, pose_shared_ptr pose_) const;




	//! builds all peptide bonds between carbon-alphas of first_res and last_res. Returns rebuild atoms vector
	POSE::atom_shared_ptr_vector build_peptide_bonds(const int first_res, const int last_res, pose_shared_ptr pose_) const;

	//! returns rebuild atoms vector
	POSE::atom_shared_ptr_vector build_fragment(const int first_res, pose_shared_ptr pose_) const;
	//! returns rebuild atoms vector
	POSE::atom_shared_ptr_vector build_backbone( pose_shared_ptr pose_) const;


};

alphabet_bb_builder_shared_ptr new_alphabet_bb_builder();


}
}
}



#endif /* ALPHABET_BB_BUILDER_H_ */
