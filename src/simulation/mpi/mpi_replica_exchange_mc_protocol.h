/*
 * mpi_replica_exchange_ca_mc_protocol.h
 *
 *  Created on: 30 Sep 2013
 *      Author: jmacdona
 */
#ifndef MPI_REPLICA_EXCHANGE_MC_PROTOCOL_H_
#define MPI_REPLICA_EXCHANGE_MC_PROTOCOL_H_

#include "simulation/mc_protocol_interface.h"
#include <limits>
#include "potentials/potentials_factory.h"
#include <map>
#include <algorithm>


#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>


namespace PRODART {
namespace POSE {
namespace SIM {
namespace MPI {





class replica_spec_element{
public:
	double beta;
	std::string potential_name;
	POTENTIALS::potentials_container_shared_ptr pot_ptr;
};

typedef std::vector<replica_spec_element> replica_spec_element_vec;
typedef std::map<int,int> int_int_map;

class exchange_request_msg{

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){
    	ar & source_replica_number;
    	ar & dest_replica_number;
    	ar & source_beta;
    	ar & source_energy;
    	ar & dest_beta;
    	ar & dest_energy;
    	ar & source_step;
    	ar & dest_step;
    	ar & accepted;
    	ar & rejected;
    	ar & confirm;
    	ar & source_rank;
    	ar & dest_rank;


    }



public:

    exchange_request_msg() : source_replica_number(999),
    dest_replica_number(999),
    source_beta(999),
    source_energy(999),
    dest_beta(999),
    dest_energy(999),
    source_step(0),
    dest_step(0),
    accepted(false),
    rejected(false),
    confirm(false),
    source_rank(-1),
    dest_rank(-1){
    }

    exchange_request_msg(int source_replica_number_,
    		int dest_replica_number_,
    		double source_beta_,
    		double source_energy_,
    		double dest_beta_,
    		double dest_energy_,
    		unsigned long source_step_,
    		int source_rank_) : source_replica_number(source_replica_number_),
    dest_replica_number(dest_replica_number_),
    source_beta(source_beta_),
    source_energy(source_energy_),
    dest_beta(dest_beta_),
    dest_energy(dest_energy_),
    source_step(source_step_),
    dest_step(0),
    accepted(false),
    rejected(false),
    confirm(false),
    source_rank(source_rank_),
    dest_rank(-1){
    }

    void flip(){
    	std::swap(source_replica_number, dest_replica_number);
    	std::swap(source_beta, dest_beta);
    	std::swap(source_energy, dest_energy);
    	std::swap(source_step, dest_step);
    	std::swap(source_rank, dest_rank);

    }

    void flip_and_accept(){
    	flip();
    	accepted = true;
    	rejected = false;
    }

    void flip_and_reject(){
    	flip();
    	accepted = false;
    	rejected = true;
    }


    void flip_and_confirm(){
    	flip();
    	confirm = true;
    }

    int source_replica_number;
    int dest_replica_number;
    double source_beta;
    double source_energy;
    double dest_beta;
    double dest_energy;
    unsigned long source_step;
    unsigned long dest_step;
    bool accepted;
    bool rejected;
    bool confirm;
    int source_rank;
    int dest_rank;

    friend std::ostream& operator<<(std::ostream& os, const exchange_request_msg& data){
        os 		<< "source_replica_number:" << data.source_replica_number << " "
        		<< "dest_replica_number:" << data.dest_replica_number << " "
        		<< "source_beta:" << data.source_beta << " "
        		<< "source_energy:" << data.source_energy << " "
        		<< "dest_beta:" << data.dest_beta << " "
        		<< "dest_energy:" << data.dest_energy << " "
        		<< "source_step:" << data.source_step << " "
        		<< "dest_step:" << data.dest_step << " "
        		<< "source_rank:" << data.source_rank << " "
        		<< "dest_rank:" << data.dest_rank << " "
        		;
        		//<< std::endl;
        if (data.accepted){
        	os << "accepted:true ";
        }
        else {
        	os << "accepted:false ";
        }

        if (data.rejected){
        	os << "rejected:true ";
        }
        else {
        	os << "rejected:false ";
        }

        if (data.confirm){
        	os << "confirm:true ";
        }
        else {
        	os << "confirm:false ";
        }

        return os;
    }
};

class mpi_replica_exchange_mc_protocol;
typedef boost::shared_ptr<mpi_replica_exchange_mc_protocol> mpi_replica_exchange_mc_protocol_shared_ptr;
mpi_replica_exchange_mc_protocol_shared_ptr new_mpi_replica_exchange_mc_protocol( boost::mpi::communicator comm,
		unsigned long steps,
		MTRand::MTRand_shared_ptr mt_ptr);

class mpi_replica_exchange_mc_protocol : public mc_protocol_interface {

	friend mpi_replica_exchange_mc_protocol_shared_ptr new_mpi_replica_exchange_mc_protocol( boost::mpi::communicator comm,
			unsigned long steps,
			MTRand::MTRand_shared_ptr mt_ptr);

private:



protected:

	replica_spec_element_vec replica_specs;
	boost::mpi::communicator rep_comm;
	int_int_map replica_rank_map;
	//int my_replica_number;
	unsigned long exchange_frequency;

	mpi_replica_exchange_mc_protocol(boost::mpi::communicator comm,
			unsigned long steps,
			double _beta,
			MTRand::MTRand_shared_ptr mt_ptr) : mc_protocol_interface(steps, _beta, mt_ptr),
			rep_comm(comm),
			//my_replica_number((int)comm.rank()),
			exchange_frequency(4000),
			exchange_msg_wait_flag(false),
			exchange_msg_confirm_flag(false),
			replica_request_energy(0),
			replica_request_step(0){

		exch_swap_fwd_requests.resize((unsigned int)rep_comm.size());
		exch_swap_result_requests.resize((unsigned int)rep_comm.size());
		rec_swap_result_vals.resize((unsigned int)rep_comm.size());


		//my_replica_number = (int)comm.rank();
		/*
		this->steps_to_run = steps;
		this->beta = _beta;
		*/
	}

	typedef boost::shared_ptr<boost::mpi::request> request_shared_ptr;
	typedef std::vector<request_shared_ptr> request_shared_ptr_vec;

	request_shared_ptr exch_swap_request;
	request_shared_ptr exch_accept_request;
	request_shared_ptr exch_confirm_request;
	request_shared_ptr_vec exch_swap_fwd_requests;
	request_shared_ptr_vec exch_swap_result_requests;

	typedef exchange_request_msg exch_receieve_type;
	typedef std::vector<exch_receieve_type> exch_receieve_type_vec;
	exch_receieve_type_vec rec_swap_result_vals;

	// don't allow exchanges while waiting
	bool exchange_msg_wait_flag;

	// don't allow exchanges while waiting
	bool exchange_msg_confirm_flag;

	POSE::atom_shared_ptr_vector3d_map replica_request_coords_backup;
	double replica_request_energy;
	unsigned long replica_request_step;



public:

	enum MessageType {ERROR=1, TEST=2, STOP_NOTIFICATION=3, SWAP_SEND=4, SWAP_RESULT=5};

	bool load_replica_specs(std::string filename);
	int get_num_replicas() const;
	unsigned long get_exchange_frequency() const {
		return exchange_frequency;
	}
	void set_exchange_frequency(unsigned long exchangeFrequency) {
		exchange_frequency = exchangeFrequency;
	}

	// virtual function. Initialise replicas
	bool stage_initialise(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	void stage_initialise_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){}
	bool stage_continue(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	bool stage_move(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	bool stage_accept_test(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	bool stage_accepted(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	void stage_accepted_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){}
	bool stage_rejected(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	void stage_rejected_print(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers){}
	bool stage_before_next_round(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);
	bool stage_finish(default_mc_state& state,
			const PRODART::POSE::pose_shared_ptr protein,
			PRODART::POSE::META::pose_meta_shared_ptr meta_data,
			PRODART::POSE::POTENTIALS::potentials_container_shared_ptr ignored,
			PRODART::POSE::MOVERS::mover_shared_ptr movers);

};






}
}
}
}


#endif /* MPI_REPLICA_EXCHANGE_MC_PROTOCOL_H_ */
