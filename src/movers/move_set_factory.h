/*
 * move_set_factory.h
 *
 *  Created on: Oct 17, 2010
 *      Author: jmacdon
 */

#ifndef MOVE_SET_FACTORY_H_
#define MOVE_SET_FACTORY_H_
#include "mover_interface.h"

#include "prodart_env/prodart_env.h"
#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include "move_set.h"
#include "ca_movers/ca_movers.h"
#include "bb_movers/bb_movers.h"
#include <vector>

namespace PRODART {
namespace POSE {
namespace MOVERS {

typedef std::vector<bool> bool_vector;


class move_set_factory {

private:
	move_set_factory(const move_set_factory&);
	static void Init();

protected:
	move_set_factory();
	virtual ~move_set_factory();

public:
	static move_set_factory* Instance();

	move_set_shared_ptr make_preset_move_set(PRODART::POSE::META::pose_meta_shared_ptr pose_meta,
			std::string preset_label,
			const MTRand::MTRand_shared_ptr rand) const;

	move_set_shared_ptr make_preset_loop_move_set_by_residue(PRODART::POSE::META::pose_meta_shared_ptr pose_meta,
			std::string preset_label,
			const MTRand::MTRand_shared_ptr rand,
			const bool_vector& loop_mask) const;



};


}
}
}





#endif /* MOVE_SET_FACTORY_H_ */
