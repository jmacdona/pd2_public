/*
 * bb_composite_ca_mover.h
 *
 *  Created on: 4 Nov 2010
 *      Author: jmacdona
 */

#ifndef BB_COMPOSITE_CA_MOVER_H_
#define BB_COMPOSITE_CA_MOVER_H_
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include <boost/shared_ptr.hpp>

namespace PRODART {
namespace POSE {
namespace MOVERS {
namespace BB {

class bb_composite_ca_mover;

typedef boost::shared_ptr<bb_composite_ca_mover> bb_composite_ca_mover_shared_ptr;

//! composite mover
class bb_composite_ca_mover : public mover_interface {

	friend bb_composite_ca_mover_shared_ptr new_bb_composite_ca_mover(const int resNum1,
			const int resNum2,
			const unsigned long high_T_steps,
			const unsigned long anneal_steps,
			const unsigned long low_T_steps,
			const double start_beta,
			const double final_beta);



private:

	bb_composite_ca_mover();
	bb_composite_ca_mover(const int resNum1,
			const int resNum2,
			const unsigned long high_T_steps,
			const unsigned long anneal_steps,
			const unsigned long low_T_steps,
			const double start_beta,
			const double final_beta);
	void init();

	int resNum1, resNum2;
	unsigned long start_beta_steps, anneal_steps, final_beta_steps;
	double start_beta, final_beta;

	std::string ca_pot_set;
	std::string bb_pot_set;

	bool ca_as_geom_pot_changed, bb_as_geom_pot_changed;
	double ca_as_geom_pot_wt, bb_as_geom_pot_wt;

public:

	bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const;

	void propagate_start_beta(const double val){
		start_beta = val;
	}
	void propagate_final_beta(const double val){
		final_beta = val;
	}
	void propagate_ca_potential_preset(const std::string str){
		ca_pot_set = str;
	}
	void propagate_bb_potential_preset(const std::string str){
		bb_pot_set = str;
	}
	void propagate_ca_potential_weight(const POTENTIALS::potentials_name& potname, const double weight){
		if (potname == POTENTIALS::potentials_name("as_geom_pot")){
			ca_as_geom_pot_changed = true;
			ca_as_geom_pot_wt = weight;
		}
	}
	void propagate_bb_potential_weight(const POTENTIALS::potentials_name& potname, const double weight){
		if (potname == POTENTIALS::potentials_name("as_geom_pot")){
			bb_as_geom_pot_changed = true;
			bb_as_geom_pot_wt = weight;
		}
	}


};



bb_composite_ca_mover_shared_ptr new_bb_composite_ca_mover(const int resNum1,
			const int resNum2,
			const unsigned long high_T_steps,
			const unsigned long anneal_steps,
			const unsigned long low_T_steps,
			const double start_beta,
			const double final_beta);

move_set_shared_ptr bb_composite_ca_mover_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const int min_seq_sep,
		const int max_seq_sep,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta,
		const PRODART::POSE::MOVERS::bool_vector& allowed_residues);

move_set_shared_ptr bb_composite_ca_mover_move_set_factory(PRODART::POSE::META::pose_meta_shared_ptr meta_data,
		const int min_seq_sep,
		const int max_seq_sep,
		const unsigned long high_T_steps,
		const unsigned long anneal_steps,
		const unsigned long low_T_steps,
		const double start_beta,
		const double final_beta);







}
}
}
}



#endif /* BB_COMPOSITE_CA_MOVER_H_ */
