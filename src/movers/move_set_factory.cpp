/*
 * move_set_factory.cpp
 *
 *  Created on: Oct 17, 2010
 *      Author: jmacdon
 */
#include "move_set_factory.h"


using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;



namespace PRODART {
namespace POSE {
namespace MOVERS {


namespace {
move_set_factory *instance = NULL;
boost::once_flag once_flag = BOOST_ONCE_INIT;
}

void move_set_factory::Init(){
	if (!instance){
		instance = new move_set_factory();
	}
}

move_set_factory* move_set_factory::Instance(){
	/*
	if (!instance){
		instance = new move_set_factory();
	}
	*/
	boost::call_once(&move_set_factory::Init, once_flag);
	return instance;
}

move_set_factory::move_set_factory(){


	cout << "move_set_factory: Initialising..." << endl;



}

move_set_factory::~move_set_factory(){

}


move_set_shared_ptr move_set_factory::make_preset_move_set(PRODART::POSE::META::pose_meta_shared_ptr pose_meta,
		std::string preset_label,
		const MTRand::MTRand_shared_ptr rand) const{
	// TODO load preset settings from database

	move_set_shared_ptr rtn_set = new_move_set();

	pose_shared_ptr pose_ = pose_meta->get_pose();

	//bool beSilent = false;
	//const int crank_max_step = 8;
	//const int crank_min_step = 2;

	if (preset_label.compare("ca_standard") == 0){
		bool_vector loop_mask(pose_->get_residue_count(), true);
		rtn_set = this->make_preset_loop_move_set_by_residue(pose_meta, "ca_standard", rand, loop_mask);
	}
	else if (preset_label.compare("ca_standard_large") == 0){

		bool_vector loop_mask(pose_->get_residue_count(), true);
		rtn_set = this->make_preset_loop_move_set_by_residue(pose_meta, "ca_standard_large", rand, loop_mask);



	}
	else if (preset_label.compare("ca_standard_small") == 0){

		bool_vector loop_mask(pose_->get_residue_count(), true);
		rtn_set = this->make_preset_loop_move_set_by_residue(pose_meta, "ca_standard_small", rand, loop_mask);

	}
	else if (preset_label.compare("bb_composite_ca") == 0){

		PRODART::POSE::MOVERS::move_set_shared_ptr bb_comp_move_set = BB::bb_composite_ca_mover_move_set_factory(pose_meta,
				3,
				7,
				10,
				200,
				100,
				0.1,
				1.2);
				//loop_mask);
		if (bb_comp_move_set->get_move_count() != 0) rtn_set->add_move(bb_comp_move_set, 1.0);

		/*
		bool_vector loop_mask(pose_->get_residue_count(), true);
		rtn_set = this->make_preset_loop_move_set_by_residue(pose_meta, "bb_composite_ca", rand, loop_mask);
		*/

	}
	else {
		cerr << "move_set_factory: move set preset no found:\t" << preset_label;
	}

	rtn_set->propagate_rand_num_gen(rand);
	return rtn_set;
}


move_set_shared_ptr move_set_factory::make_preset_loop_move_set_by_residue(PRODART::POSE::META::pose_meta_shared_ptr pose_meta,
		std::string preset_label,
		const MTRand::MTRand_shared_ptr rand,
		const bool_vector& loop_mask) const{
	// TODO load preset settings from database

	move_set_shared_ptr rtn_set = new_move_set();

	pose_shared_ptr pose_ = pose_meta->get_pose();

	bool beSilent = true;
#ifndef NDEBUG
	beSilent = false;
#endif
	const int crank_max_step = 8;
	const int crank_min_step = 2;

	if (preset_label.compare("ca_standard") == 0){
		rtn_set->add_move(this->make_preset_loop_move_set_by_residue(pose_meta, "ca_standard_small", rand, loop_mask), 1.0);
		rtn_set->add_move(this->make_preset_loop_move_set_by_residue(pose_meta, "ca_standard_large", rand, loop_mask), 0.1);
	}
	else if (preset_label.compare("ca_standard_large") == 0){

		PRODART::POSE::MOVERS::move_set_shared_ptr dih_move_set = PRODART::POSE::MOVERS::CA::ca_single_dihedral_uni_dist_move_set_factory(pose_meta, PI / 2.0,
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr ang_move_set = PRODART::POSE::MOVERS::CA::ca_single_angle_uni_dist_move_set_factory(pose_meta, 10.0 * (PI / 180.0),
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr pinch_move_set = PRODART::POSE::MOVERS::CA::ca_local_angle_pinch_uni_dist_move_set_factory(pose_meta, 20.0 * (PI / 180.0),
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr bond_move_set = PRODART::POSE::MOVERS::CA::ca_local_bond_uni_dist_move_set_factory(pose_meta, 1.0,
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr crank_move_set = PRODART::POSE::MOVERS::CA::ca_single_crankshaft_uni_dist_move_set_factory(pose_meta,
				PI / 2.0,
				crank_min_step,
				crank_max_step,
				loop_mask);

		if (dih_move_set->get_move_count() != 0) rtn_set->add_move(dih_move_set, 1.0);
		if (ang_move_set->get_move_count() != 0) rtn_set->add_move(ang_move_set, 1.0);
		if (bond_move_set->get_move_count() != 0) rtn_set->add_move(bond_move_set, 1.0);
		if (crank_move_set->get_move_count() != 0) rtn_set->add_move(crank_move_set, 2.0);
		if (pinch_move_set->get_move_count() != 0) rtn_set->add_move(pinch_move_set, 1.0);

		if (!beSilent){

			cout << "dih movers made: " << dih_move_set->get_move_count() << endl;
			cout << "ang movers made: " << ang_move_set->get_move_count() << endl;
			cout << "bond movers made: " << bond_move_set->get_move_count() << endl;
			cout << "crank movers made: " << crank_move_set->get_move_count() << endl;
			cout << "pinch movers made: " << pinch_move_set->get_move_count() << endl;

		}

	}
	else if (preset_label.compare("ca_standard_small") == 0){


		//note: change this to get resCount of longest chain
		int resCount = 0;//protein.getResidueCount();
		int long_chain = 0;
		const int chain_count = pose_->get_chain_count();
		for (int i = 0; i < chain_count; i++){
			if (pose_->get_chain(i)->length() > resCount ){
				resCount = pose_->get_chain(i)->length();
				long_chain = i;
			}
		}

		const int half_resCount = 1 + (resCount / 2);

		const double length_per_res = 3.8;

		const double move_margin = ENV::get_option_value<double>("sim:ca:move_margin") ;//1.5;
		const unsigned long nb_pair_update_freq = 30;

		const double max_cartesian_move = move_margin / static_cast<double>(nb_pair_update_freq);

		// from random walk - ideal polymer
		const  double expected_length = length_per_res * (2.0 + std::sqrt(static_cast<double>(half_resCount)));

		double real_length = (pose_->get_chain(long_chain)->get_ca_pos(0)
				- pose_->get_chain(long_chain)->get_ca_pos(half_resCount)).mod();

		if (!beSilent){
			cout << "\nexpected_length:\t" << expected_length << endl;
			cout << "real_length\t" << real_length << endl;
		}
		const double max_angle = std::asin(max_cartesian_move
			/ expected_length);
		//added 7.5 / 5.0 multiplying factor
		const  double expected_crank_length = length_per_res * (2.0 + std::sqrt(static_cast<double>(crank_max_step/2)));

		const double max_crank_angle = std::asin(max_cartesian_move
			/ expected_crank_length);

		PRODART::POSE::MOVERS::move_set_shared_ptr dih_move_set = PRODART::POSE::MOVERS::CA::ca_single_dihedral_uni_dist_move_set_factory(pose_meta, max_angle,
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr ang_move_set = PRODART::POSE::MOVERS::CA::ca_single_angle_uni_dist_move_set_factory(pose_meta, max_angle,
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr pinch_move_set = PRODART::POSE::MOVERS::CA::ca_local_angle_pinch_uni_dist_move_set_factory(pose_meta, max_angle,
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr bond_move_set = PRODART::POSE::MOVERS::CA::ca_local_bond_uni_dist_move_set_factory(pose_meta, max_cartesian_move,
				loop_mask);
		PRODART::POSE::MOVERS::move_set_shared_ptr crank_move_set = PRODART::POSE::MOVERS::CA::ca_single_crankshaft_uni_dist_move_set_factory(pose_meta,
				max_crank_angle,
				crank_min_step,
				crank_max_step,
				loop_mask);

		if (dih_move_set->get_move_count() != 0) rtn_set->add_move(dih_move_set, 1.0);
		if (ang_move_set->get_move_count() != 0) rtn_set->add_move(ang_move_set, 1.0);
		if (bond_move_set->get_move_count() != 0) rtn_set->add_move(bond_move_set, 0.5);
		if (crank_move_set->get_move_count() != 0) rtn_set->add_move(crank_move_set, 2.0);
		if (pinch_move_set->get_move_count() != 0) rtn_set->add_move(pinch_move_set, 1.0);



		if (!beSilent){

			cout << "dih movers made: " << dih_move_set->get_move_count() << endl;
			cout << "ang movers made: " << ang_move_set->get_move_count() << endl;
			cout << "bond movers made: " << bond_move_set->get_move_count() << endl;
			cout << "crank movers made: " << crank_move_set->get_move_count() << endl;
			cout << "pinch movers made: " << pinch_move_set->get_move_count() << endl;

			double crankMoveCutoff = max_crank_angle;
			double dihedralMoveCutoff = max_angle;
			double angleMoveCutoff = max_angle;
			double cartesianMoveCutoff = max_cartesian_move;
			cout << endl
				<< "crankMoveCutoff:\t" << crankMoveCutoff * (180.0/PI) << endl
				<< "dihedralMoveCutoff\t" << dihedralMoveCutoff * (180.0/PI) << endl
				<< "angleMoveCutoff\t" << angleMoveCutoff * (180.0/PI) << endl
				<< "cartesianMoveCutoff\t" << cartesianMoveCutoff << endl
				<< endl;
		}


	}
	else if (preset_label.compare("bb_composite_ca") == 0){
		PRODART::POSE::MOVERS::move_set_shared_ptr bb_comp_move_set = BB::bb_composite_ca_mover_move_set_factory(pose_meta,
				1,
				8,
				10,
				200,
				100,
				0.1,
				1.2,
				loop_mask);
		if (bb_comp_move_set->get_move_count() != 0) rtn_set->add_move(bb_comp_move_set, 1.0);
	}
	else {
		cerr << "move_set_factory: move set preset no found:\t" << preset_label;
	}

	rtn_set->propagate_rand_num_gen(rand);
	return rtn_set;
}



}
}
}

