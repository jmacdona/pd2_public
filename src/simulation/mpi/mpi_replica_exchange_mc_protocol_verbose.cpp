/*
 * mpi_replica_exchange_mc_protocol_verbose.cpp
 *
 *  Created on: 3 Oct 2013
 *      Author: jmacdona
 */
#include "mpi_replica_exchange_mc_protocol_verbose.h"




namespace PRODART {
namespace POSE {
namespace SIM {
namespace MPI {

mpi_replica_exchange_mc_protocol_verbose_shared_ptr new_mpi_replica_exchange_mc_protocol_verbose( boost::mpi::communicator comm,
		unsigned long steps,
		MTRand::MTRand_shared_ptr mt_ptr){
	mpi_replica_exchange_mc_protocol_verbose_shared_ptr ptr(new mpi_replica_exchange_mc_protocol_verbose(comm, steps, 1.0, mt_ptr));
	return ptr;
}



}
}
}
}
