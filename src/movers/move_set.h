//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * move_set.h
 *
 *  Created on: Mar 1, 2010
 *      Author: jmacdon
 */

#ifndef MOVE_SET_H_
#define MOVE_SET_H_
#include "mover_interface.h"

#include "prodart_env/prodart_env.h"
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


class move_set_entry;
class move_set;

typedef boost::shared_ptr<move_set> move_set_shared_ptr;
typedef std::vector<move_set_entry> move_set_entry_vector;

class move_set : public mover_interface {

	friend move_set_shared_ptr new_move_set();

public:

	virtual ~move_set(){
	}

	bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data	) const;

	void add_move(mover_shared_ptr ptr, const double weight);

	void calc_prob_accepts();

	void propagate_rand_num_gen(const MTRand::MTRand_shared_ptr ptr);

	void propagate_start_beta(const double val);
	void propagate_final_beta(const double val);
	void propagate_ca_potential_preset(const std::string str);
	void propagate_bb_potential_preset(const std::string str);
	void propagate_ca_potential_weight(const POTENTIALS::potentials_name&, const double weight);
	void propagate_bb_potential_weight(const POTENTIALS::potentials_name&, const double weight);

	unsigned int get_move_count() const;

private:

	move_set();

	double find_max_weight() const;

	move_set_entry_vector move_vec;

	//MTRand::MTRand_shared_ptr rand_gen;

	static int max_tries;


};

move_set_shared_ptr new_move_set();


class move_set_entry{
public:
	mover_shared_ptr mover;
	double weight;
	double prob_accept;
};

class move_set_bad_move_add: public std::exception{
  virtual const char* what() const throw()
  {
    return "move_set_bad_move_add";
  }
};



}
}
}






#endif /* MOVE_SET_H_ */
