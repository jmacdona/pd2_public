/*
 * ca_radgyr.cpp
 *
 *  Created on: 8 Sep 2010
 *      Author: jmacdona
 */

#include "ca_radgyr.h"

using namespace boost;
using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE::POTENTIALS;
using namespace PRODART::POSE::META;
using namespace PRODART::POSE_UTILS;
using namespace std;


namespace PRODART {
namespace POSE {
namespace POTENTIALS {
namespace CA{

potential_shared_ptr new_ca_radgyr(){
	potential_shared_ptr ptr(new ca_radgyr());
	return ptr;
}

double ca_radgyr::mean_native_radgyr(const double resCount) const{
	return native_mean_a0 + native_mean_a1 * pow(resCount,(2.0/5.0) );
}
double ca_radgyr::mean_random_walk_radgyr(const double resCount) const{
	return rand_walk_mean_a0 + rand_walk_mean_a1 * pow(resCount,(3.0/5.0) );
}

double ca_radgyr::stddev_native_radgyr(const double resCount) const{
	return native_stddev_a0 + native_stddev_a1 * pow(resCount,(2.0/5.0) );
}
double ca_radgyr::stddev_random_walk_radgyr(const double resCount) const{
	return rand_walk_stddev_a0 + rand_walk_stddev_a1 * pow(resCount,(3.0/5.0) );
}


inline double ca_radgyr::k_native_radgyr(const double resCount) const{
	return 1.0 / (2.0 * stddev_native_radgyr(resCount));
}

inline double ca_radgyr::k_random_walk_radgyr(const double resCount) const{
	return 1.0 / (2.0 * stddev_random_walk_radgyr(resCount));
}
double ca_radgyr::mean_native_radgyr( const PRODART::POSE::const_pose_shared_ptr protein, const int chain) const{
	const double numRes = static_cast<double>(protein->get_chain(chain)->length());
	return mean_native_radgyr(numRes);
}

double ca_radgyr::stddev_native_radgyr( const PRODART::POSE::const_pose_shared_ptr protein, const int chain) const{
	const double numRes = static_cast<double>(protein->get_chain(chain)->length());
	return stddev_native_radgyr(numRes);
}


double ca_radgyr::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	double total_energy = 0;
	const const_pose_shared_ptr protein = _pose_meta->get_pose();
	const ca_pose_meta_shared_ptr ca_meta_dat = static_pointer_cast<ca_pose_meta, pose_meta_interface>(_pose_meta);
	const int numChains = protein->get_chain_count();

	for (int i = 0; i < numChains; i++){
		const_chain_shared_ptr currChain = protein->get_chain(i);
		const double numRes = static_cast<double>(currChain->length());

		const double k_nat = k_native_radgyr(numRes);
		const double l0_nat = mean_native_radgyr(numRes);

		const double k_rw = k_random_walk_radgyr(numRes);
		const double l0_rw = mean_random_walk_radgyr(numRes);

		const double radgyr = PRODART::POSE_UTILS::get_ca_radgyr(protein, i);//PRODART::PROT::CA_radgyr(protein, i);

		total_energy += k_nat * pow((radgyr - l0_nat),2)
						- k_rw * pow((radgyr - l0_rw),2);

	}

	return energies_map.add_energy_component(this->name_vector[0], total_energy);
}

double ca_radgyr::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr _pose_meta,
		potentials_energies_map& energies_map) const{
	return get_energy(_pose_meta, energies_map);
}

bool ca_radgyr::init(){

	return true;
}

//! does it provide energies labeled with this name
/*
bool ca_radgyr::provides(const potentials_name& query_name) const{
	if (query_name == name){
		return true;
	}

	return false;
}
*/

/*
std::istream& ca_radgyr::load_data( std::istream& input ){

}
*/

ca_radgyr::ca_radgyr() {

	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("ca_radgyr"));

	native_mean_a0 = 1.51669;
	native_mean_a1 = 1.76209;

	rand_walk_mean_a0 = -2.64485;//-10.8;//-4.8;//-2.2 ;//-1.7;//			orig: -1.5337;
	rand_walk_mean_a1 = 1.91563;//2.5;// 2.0;//1.8 ;//1.7;//			orig:  1.58684;

	native_stddev_a0 = 0.21171;
	native_stddev_a1 = 0.0236863;

	rand_walk_stddev_a0 = -0.20;		//    			orig: -0.00223403;
	rand_walk_stddev_a1 = 0.143831;//	orig:  0.123831;
}




}
}
}
}

