//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * mover_defs.h
 *
 *  Created on: Jun 6, 2010
 *      Author: jmacdon
 */

#ifndef MOVER_DEFS_H_
#define MOVER_DEFS_H_

namespace PRODART {
namespace POSE {
namespace MOVERS {

class mover_flags {

public:
	bool move_completed;
	bool is_large_move;
	//! largest distance moved by any atom
	double move_dist;

public:
	mover_flags(){
		move_completed = true;
		is_large_move = true;
		move_dist = 0;
	}


};

}
}
}

#endif /* MOVER_DEFS_H_ */
