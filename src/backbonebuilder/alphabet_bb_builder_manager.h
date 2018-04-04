/*
 * alphabet_bb_builder_manager.h
 *
 *  Created on: Jun 7, 2012
 *      Author: jmacdona
 */

#ifndef ALPHABET_BB_BUILDER_MANAGER_H_
#define ALPHABET_BB_BUILDER_MANAGER_H_

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

#include "alphabet_bb_builder.h"

namespace PRODART {
namespace POSE {
namespace BB_BUILDER {

class alphabet_bb_builder_manager;

typedef boost::shared_ptr<alphabet_bb_builder_manager> alphabet_bb_builder_manager_shared_ptr;
typedef boost::shared_ptr<const alphabet_bb_builder_manager> const_alphabet_bb_builder_manager_shared_ptr;

typedef std::map<int, const_alphabet_bb_builder_shared_ptr> int_const_alphabet_bb_builder_shared_ptr_map;


class alphabet_bb_builder_manager {



private:

	alphabet_bb_builder_manager(const alphabet_bb_builder_manager&);

	static void Init();

protected:
	alphabet_bb_builder_manager();
	//virtual ~alphabet_bb_builder_manager();

	int_const_alphabet_bb_builder_shared_ptr_map alphabets;

	int default_even_len,
		default_odd_len,
		default_len;

public:

	static const_alphabet_bb_builder_manager_shared_ptr Instance();

	virtual ~alphabet_bb_builder_manager();

	const_alphabet_bb_builder_shared_ptr get(const int alphabet_len) const;

	POSE::atom_shared_ptr_vector build_backbone( pose_shared_ptr pose_) const;

	//! builds all peptide bonds between carbon-alphas of first_res and last_res. Returns rebuild atoms vector
	POSE::atom_shared_ptr_vector build_peptide_bonds(const int first_res, const int last_res, pose_shared_ptr pose_) const;


};






}
}
}



#endif /* ALPHABET_BB_BUILDER_MANAGER_H_ */
