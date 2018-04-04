/*
 * strand_pairing.cpp
 *
 *  Created on: 28 Jan 2011
 *      Author: jmacdona
 */

#include "strand_pairing_rst.h"

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

const double ca_lower_cutoff = 4.4;
const double ca_upper_cutoff = 5.5;
const double ca_next_lower_cutoff = 5.4,
		ca_next_upper_cutoff = 7.5;
const double ca_next_weight = 0.5;

const double cb_lower_cutoff = 3.8;
const double cb_upper_cutoff = 6.2;

const double bb_lower_cutoff = 1.6;
const double bb_upper_cutoff = 2.2;
const double bb_angle_cutoff = degrees_to_radians(90.0);  //originally set to 60 degrees not sure why.

const double penalty = 500.0;

potential_shared_ptr new_strand_pairing_rst(){
	potential_shared_ptr ptr(new strand_pairing_rst());
	return ptr;
}

strand_pairing_rst::strand_pairing_rst() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("strand_pairing_rst"));
}



bool strand_pairing_rst::init(){

	return true;
}

double strand_pairing_rst::get_dist_energy(const double lower_cutoff,
		const double upper_cutoff,
		const double dist) const{


	if (dist < lower_cutoff ){
		return std::pow(dist - lower_cutoff, 2);
	}
	else if (dist > upper_cutoff){
		return std::pow(dist - upper_cutoff, 2);
	}
	else {
		return 0;
	}

}

double strand_pairing_rst::get_dE_ddist(const double lower_cutoff,
		const double upper_cutoff,
		const double dist) const{
	if (dist < lower_cutoff ){
		return 2.0 * (dist - lower_cutoff);
	}
	else if (dist > upper_cutoff){
		return 2.0 * (dist - upper_cutoff);
	}
	else {
		return 0;
	}
}


double strand_pairing_rst::get_pair_energy(const double lower_cutoff,
		const double upper_cutoff, POSE::atom_shared_ptr at1, POSE::atom_shared_ptr at2,
		const double weight) const{
	const double dist = (at1->get_coords() - at2->get_coords()).mod();
	return weight * get_dist_energy(lower_cutoff, upper_cutoff, dist);
}

double strand_pairing_rst::get_pair_energy_grad(const double lower_cutoff,
		const double upper_cutoff, POSE::atom_shared_ptr at1, POSE::atom_shared_ptr at2,
		const double weight,
		PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const{
	const int seq_num1 = at1->get_seq_num();
	const int seq_num2 = at2->get_seq_num();
	const double dist = (at1->get_coords() - at2->get_coords()).mod();
	const vector3d grad_vec = (at1->get_coords() - at2->get_coords()) / dist;
	const double dE_dd = get_dE_ddist(lower_cutoff, upper_cutoff, dist);
	grad[seq_num1] += weight * pot_weight * dE_dd * grad_vec;
	grad[seq_num2] += -weight * pot_weight * dE_dd * grad_vec;
	return weight * get_dist_energy(lower_cutoff, upper_cutoff, dist);
}

double strand_pairing_rst::get_ca_pair_energy_grad(POSE::pose_shared_ptr pose_,
		const int res1, const int res2,
		const bool is_parallel, const double weight,
		PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const{
	atom_shared_ptr at1 = pose_->get_bb_atom(POSE::CA, res1);
	atom_shared_ptr at2 = pose_->get_bb_atom(POSE::CA, res2);
	if (at1 && at2){
		if (at1->isActiveAndSet() && at2->isActiveAndSet()){
			double energy =  get_pair_energy_grad(ca_lower_cutoff, ca_upper_cutoff, at1, at2, weight, grad, pot_weight);

			residue_shared_ptr res1_ptr = pose_->get_residue(res1);
			residue_shared_ptr res2_ptr = pose_->get_residue(res2);

			residue_shared_ptr res1_prev_next[2] = {res1_ptr->get_prev_residue(), res1_ptr->get_next_residue()};
			residue_shared_ptr res2_prev_next[2] = {res2_ptr->get_prev_residue(), res2_ptr->get_next_residue()};

			for (int i = 0; i < 2; i++){
				if (res2_prev_next[i]){
					if (res2_prev_next[i]->get_chain() == res2_ptr->get_chain()){
						atom_shared_ptr this_ca = res2_prev_next[i]->get_bb_atom(POSE::CA);
						if (this_ca->isActiveAndSet() && this_ca != at1){
							energy += this->get_pair_energy_grad(ca_next_lower_cutoff, ca_next_upper_cutoff, at1, this_ca, ca_next_weight * weight, grad, pot_weight);
						}
					}
				}
				if (res1_prev_next[i]){
					if (res1_prev_next[i]->get_chain() == res1_ptr->get_chain()){
						atom_shared_ptr this_ca = res1_prev_next[i]->get_bb_atom(POSE::CA);
						if (this_ca->isActiveAndSet() && this_ca != at2){
							energy += this->get_pair_energy_grad(ca_next_lower_cutoff, ca_next_upper_cutoff, at2, this_ca, ca_next_weight * weight, grad, pot_weight);
						}
					}
				}
			}
			//assert(is_vector3d_vector_valid(grad));
			return energy;

			/*
			const int seq_num1 = at1->get_seq_num();
			const int seq_num2 = at2->get_seq_num();
			const double dist = (at1->get_coords() - at2->get_coords()).mod();
			const vector3d grad_vec = (at1->get_coords() - at2->get_coords()) / dist;
			const double dE_dd = get_dE_ddist(ca_lower_cutoff, ca_upper_cutoff, dist);
			grad[seq_num1] += weight * pot_weight * dE_dd * grad_vec;
			grad[seq_num2] += -weight * pot_weight * dE_dd * grad_vec;
			return weight * get_dist_energy(ca_lower_cutoff, ca_upper_cutoff, dist);
			*/
		}
	}
	return 0;
}
double strand_pairing_rst::get_cb_pair_energy_grad(POSE::pose_shared_ptr pose_,
		const int res1, const int res2,
		const bool is_parallel, const double weight,
		PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const{
	atom_shared_ptr at1 = pose_->get_bb_atom(POSE::CB, res1);
	atom_shared_ptr at2 = pose_->get_bb_atom(POSE::CB, res2);
	if (at1 && at2){
		if (at1->isActiveAndSet() && at2->isActiveAndSet()){
			const double energy =  get_pair_energy_grad(cb_lower_cutoff, cb_upper_cutoff, at1, at2, weight, grad, pot_weight);


			return energy;

		}
	}
	return 0;
}

double strand_pairing_rst::get_bb_pair_energy_grad(POSE::pose_shared_ptr pose_,
		const int res1, const int res2,
		const bool is_parallel, const double weight,
		PRODART::UTILS::vector3d_vector& grad, const double pot_weight) const{
	if (!is_parallel){
		atom_shared_ptr at1_h = pose_->get_bb_atom(POSE::H, res1);
		atom_shared_ptr at2_h = pose_->get_bb_atom(POSE::H, res2);
		atom_shared_ptr at1_n = pose_->get_bb_atom(POSE::N, res1);
		atom_shared_ptr at2_n = pose_->get_bb_atom(POSE::N, res2);
		atom_shared_ptr at1_o = pose_->get_bb_atom(POSE::O, res1);
		atom_shared_ptr at2_o = pose_->get_bb_atom(POSE::O, res2);
		atom_shared_ptr at1_c = pose_->get_bb_atom(POSE::C, res1);
		atom_shared_ptr at2_c = pose_->get_bb_atom(POSE::C, res2);
		//cout << "db0\n";
		if (at1_h && at2_h
				&& at1_n && at2_n
				&& at1_o && at2_o
				&& at1_c && at2_c){
			//cout << "db1\n";
			if (at1_h->isActiveAndSet() && at2_h->isActiveAndSet()
					&& at1_n->isActiveAndSet() && at2_n->isActiveAndSet()
					&& at1_o->isActiveAndSet() && at2_o->isActiveAndSet()
					&& at1_c->isActiveAndSet() && at2_c->isActiveAndSet()){
				//cout << "db2\n";

				//cout << "db3\n";
				//vector3d n1_h1_vec = at1_h->get_coords() - at1_n->get_coords();
				//vector3d n2_h2_vec = at2_h->get_coords() - at2_n->get_coords();

				//vector3d c1_o1_vec = at1_o->get_coords() - at1_c->get_coords();
				//vector3d c2_o2_vec = at2_o->get_coords() - at2_c->get_coords();

				//if ((n1_h1_vec.dot(c2_o2_vec) < 0) && (c1_o1_vec.dot(n2_h2_vec) < 0)){
				if (angle(at1_n->get_coords(), at1_h->get_coords(), at2_o->get_coords()) > bb_angle_cutoff
						&& angle(at1_h->get_coords(), at2_o->get_coords(), at2_c->get_coords()) > bb_angle_cutoff
						&& angle(at2_n->get_coords(), at2_h->get_coords(), at1_o->get_coords()) > bb_angle_cutoff
						&& angle(at2_h->get_coords(), at1_o->get_coords(), at1_c->get_coords()) > bb_angle_cutoff){
					//cout << "db4\n";
					//|| (n2_h2_vec.dot(c1_o1_vec) < 0 && c2_o2_vec.dot(n1_h1_vec) < 0) ){
					//const double h1_o2_dist = (at1_h->get_coords() - at2_o->get_coords()).mod();
					//const double h2_o1_dist = (at2_h->get_coords() - at1_o->get_coords()).mod();

					const double enrg = get_pair_energy_grad(bb_lower_cutoff, bb_upper_cutoff, at1_h, at2_o, weight, grad, pot_weight)
									+ get_pair_energy_grad(bb_lower_cutoff, bb_upper_cutoff, at2_h, at1_o, weight, grad, pot_weight);
					/*
					{
						const int seq_num1 = at1_h->get_seq_num();
						const int seq_num2 = at2_o->get_seq_num();
						const double dist = h1_o2_dist;//(at1->get_coords() - at2->get_coords()).mod();
						const vector3d grad_vec = (at1_h->get_coords() - at2_o->get_coords()) / dist;
						const double dE_dd = get_dE_ddist(bb_lower_cutoff, bb_upper_cutoff, dist);
						grad[seq_num1] += weight * pot_weight * dE_dd * grad_vec;
						grad[seq_num2] += -weight * pot_weight * dE_dd * grad_vec;
					}
					{
						const int seq_num1 = at2_h->get_seq_num();
						const int seq_num2 = at1_o->get_seq_num();
						const double dist = h2_o1_dist;//(at1->get_coords() - at2->get_coords()).mod();
						const vector3d grad_vec = (at2_h->get_coords() - at1_o->get_coords()) / dist;
						const double dE_dd = get_dE_ddist(bb_lower_cutoff, bb_upper_cutoff, dist);
						grad[seq_num1] += weight * pot_weight * dE_dd * grad_vec;
						grad[seq_num2] += -weight * pot_weight * dE_dd * grad_vec;
					}



					const double enrg = (weight * get_dist_energy(bb_lower_cutoff, bb_upper_cutoff, h1_o2_dist))
							+ (weight * get_dist_energy(bb_lower_cutoff, bb_upper_cutoff, h2_o1_dist)) ;
					 */
					//assert(is_vector3d_vector_valid(grad));
					return enrg;
				}
				else {
					//cout << "db5\n";
					return penalty; // + this->get_ca_pair_energy_grad(pose_, res1, res2, is_parallel, 1.0 * weight, grad, pot_weight);
				}


			} //
		}
	}
	else {
		//PARALLEL STRANDS
		atom_shared_ptr at1_h = pose_->get_bb_atom(POSE::H, res1);
		atom_shared_ptr at2_h = pose_->get_bb_atom(POSE::H, res2);
		atom_shared_ptr at1_n = pose_->get_bb_atom(POSE::N, res1);
		atom_shared_ptr at2_n = pose_->get_bb_atom(POSE::N, res2);
		atom_shared_ptr at1_o = pose_->get_bb_atom(POSE::O, res1);
		atom_shared_ptr at2_o = pose_->get_bb_atom(POSE::O, res2);
		atom_shared_ptr at1_c = pose_->get_bb_atom(POSE::C, res1);
		atom_shared_ptr at2_c = pose_->get_bb_atom(POSE::C, res2);

		const const_residue_shared_ptr res1_ptr = pose_->get_residue(res1);
		const const_residue_shared_ptr res2_ptr = pose_->get_residue(res2);

		const const_residue_shared_ptr res1_m1_ptr = res1_ptr->get_prev_residue();
		//bool res1_m1OK = false;
		atom_shared_ptr at1_m1_o = atom_shared_ptr(), at1_m1_c = atom_shared_ptr();
		if (res1_m1_ptr){
				//res1_m1OK = true;
				at1_m1_o = pose_->get_bb_atom(POSE::O, res1-1);
				at1_m1_c = pose_->get_bb_atom(POSE::C, res1-1);
		}

		const const_residue_shared_ptr res1_p1_ptr = res1_ptr->get_next_residue();
		//bool res1_p1OK = false;
		atom_shared_ptr at1_p1_h = atom_shared_ptr(), at1_p1_n = atom_shared_ptr();
		if (res1_p1_ptr){
				//res1_p1OK = true;
				at1_p1_h = pose_->get_bb_atom(POSE::H, res1+1);
				at1_p1_n = pose_->get_bb_atom(POSE::N, res1+1);
		}

		const const_residue_shared_ptr res2_m1_ptr = res2_ptr->get_prev_residue();
		//bool res2_m1OK = false;
		atom_shared_ptr at2_m1_o = atom_shared_ptr(), at2_m1_c = atom_shared_ptr();
		if (res2_m1_ptr){
				//res2_m1OK = true;
				at2_m1_o = pose_->get_bb_atom(POSE::O, res2-1);
				at2_m1_c = pose_->get_bb_atom(POSE::C, res2-1);
		}

		const const_residue_shared_ptr res2_p1_ptr = res2_ptr->get_next_residue();
		//bool res2_p1OK = false;
		atom_shared_ptr at2_p1_h = atom_shared_ptr(), at2_p1_n = atom_shared_ptr();
		if (res2_p1_ptr){
				//res2_p1OK = true;
				at2_p1_h = pose_->get_bb_atom(POSE::H, res2+1);
				at2_p1_n = pose_->get_bb_atom(POSE::N, res2+1);
		}

		/*
		 * Check for at1_h -> at2_m1_o && at1_o -> at2_p1_h H-bonding
		 */
		if (at1_h->isActiveAndSet() && at2_m1_o->isActiveAndSet() && at1_o->isActiveAndSet() && at2_p1_h->isActiveAndSet()
				&& at2_m1_c->isActiveAndSet() && at1_n->isActiveAndSet() && at2_p1_n->isActiveAndSet() && at1_c->isActiveAndSet()
				&& angle(at1_h->get_coords(), at2_m1_o->get_coords(), at2_m1_c->get_coords()) > bb_angle_cutoff
				&& angle(at1_n->get_coords(), at1_h->get_coords(), at2_m1_o->get_coords()) > bb_angle_cutoff
				&& angle(at1_o->get_coords(), at2_p1_h->get_coords(), at2_p1_n->get_coords()) > bb_angle_cutoff
				&& angle(at1_c->get_coords(), at1_o->get_coords(), at2_p1_h->get_coords()) > bb_angle_cutoff){
			return this->get_pair_energy_grad(bb_lower_cutoff, bb_upper_cutoff, at1_h, at2_m1_o, weight, grad, pot_weight)
					+ this->get_pair_energy_grad(bb_lower_cutoff, bb_upper_cutoff, at1_o, at2_p1_h, weight, grad, pot_weight);
		}
		/*
		 * Check for at2_h -> at1_m1_o && at2_o -> at1_p1_h H-bonding
		 */
		else if (at2_h->isActiveAndSet() && at1_m1_o->isActiveAndSet() && at2_o->isActiveAndSet() && at1_p1_h->isActiveAndSet()
				&& at1_m1_c->isActiveAndSet() && at2_n->isActiveAndSet() && at1_p1_n->isActiveAndSet() && at2_c->isActiveAndSet()
				&& angle(at2_h->get_coords(), at1_m1_o->get_coords(), at1_m1_c->get_coords()) > bb_angle_cutoff
				&& angle(at2_n->get_coords(), at2_h->get_coords(), at1_m1_o->get_coords()) > bb_angle_cutoff
				&& angle(at2_o->get_coords(), at1_p1_h->get_coords(), at1_p1_n->get_coords()) > bb_angle_cutoff
				&& angle(at2_c->get_coords(), at2_o->get_coords(), at1_p1_h->get_coords()) > bb_angle_cutoff){
			return this->get_pair_energy_grad(bb_lower_cutoff, bb_upper_cutoff, at2_h, at1_m1_o, weight, grad, pot_weight)
					+ this->get_pair_energy_grad(bb_lower_cutoff, bb_upper_cutoff, at2_o, at1_p1_h, weight, grad, pot_weight);
		}
		else {
			return penalty; // + this->get_ca_pair_energy_grad(pose_, res1, res2, is_parallel, 1.0 * weight, grad, pot_weight);
		}

	}
	return 0;
}

double strand_pairing_rst::get_ca_pair_energy(POSE::pose_shared_ptr pose_, const int res1, const int res2, const bool is_parallel, const double weight) const{
	atom_shared_ptr at1 = pose_->get_bb_atom(POSE::CA, res1);
	atom_shared_ptr at2 = pose_->get_bb_atom(POSE::CA, res2);
	if (at1 && at2){
		if (at1->isActiveAndSet() && at2->isActiveAndSet()){
			/*
			const double dist = (at1->get_coords() - at2->get_coords()).mod();
			return weight * get_dist_energy(ca_lower_cutoff, ca_upper_cutoff, dist);
			*/
			double energy =  this->get_pair_energy(ca_lower_cutoff, ca_upper_cutoff, at1, at2, weight);

			residue_shared_ptr res1_ptr = pose_->get_residue(res1);
			residue_shared_ptr res2_ptr = pose_->get_residue(res2);

			residue_shared_ptr res1_prev_next[2] = {res1_ptr->get_prev_residue(), res1_ptr->get_next_residue()};
			residue_shared_ptr res2_prev_next[2] = {res2_ptr->get_prev_residue(), res2_ptr->get_next_residue()};

			for (int i = 0; i < 2; i++){
				if (res2_prev_next[i]){
					if (res2_prev_next[i]->get_chain() == res2_ptr->get_chain()){
						atom_shared_ptr this_ca = res2_prev_next[i]->get_bb_atom(POSE::CA);
						if (this_ca->isActiveAndSet()){
							energy += this->get_pair_energy(ca_next_lower_cutoff, ca_next_upper_cutoff, at1, this_ca, ca_next_weight * weight);
						}
					}
				}
				if (res1_prev_next[i]){
					if (res1_prev_next[i]->get_chain() == res1_ptr->get_chain()){
						atom_shared_ptr this_ca = res1_prev_next[i]->get_bb_atom(POSE::CA);
						if (this_ca->isActiveAndSet()){
							energy += this->get_pair_energy(ca_next_lower_cutoff, ca_next_upper_cutoff, at2, this_ca, ca_next_weight * weight);
						}
					}
				}
			}

			return energy;
		}
	}
	return 0;
}
double strand_pairing_rst::get_cb_pair_energy(POSE::pose_shared_ptr pose_, const int res1, const int res2, const bool is_parallel, const double weight) const{
	atom_shared_ptr at1 = pose_->get_bb_atom(POSE::CB, res1);
	atom_shared_ptr at2 = pose_->get_bb_atom(POSE::CB, res2);
	if (at1 && at2){
		if (at1->isActiveAndSet() && at2->isActiveAndSet()){
			/*
			const double dist = (at1->get_coords() - at2->get_coords()).mod();
			return weight * get_dist_energy(ca_lower_cutoff, ca_upper_cutoff, dist);
			*/
			const double energy =  this->get_pair_energy(cb_lower_cutoff, cb_upper_cutoff, at1, at2, weight);

			return energy;
		}
	}
	return 0;
}
double strand_pairing_rst::get_bb_pair_energy(POSE::pose_shared_ptr pose_, const int res1, const int res2, const bool is_parallel, const double weight) const{
	if (!is_parallel){
		atom_shared_ptr at1_h = pose_->get_bb_atom(POSE::H, res1);
		atom_shared_ptr at2_h = pose_->get_bb_atom(POSE::H, res2);
		atom_shared_ptr at1_n = pose_->get_bb_atom(POSE::N, res1);
		atom_shared_ptr at2_n = pose_->get_bb_atom(POSE::N, res2);
		atom_shared_ptr at1_o = pose_->get_bb_atom(POSE::O, res1);
		atom_shared_ptr at2_o = pose_->get_bb_atom(POSE::O, res2);
		atom_shared_ptr at1_c = pose_->get_bb_atom(POSE::C, res1);
		atom_shared_ptr at2_c = pose_->get_bb_atom(POSE::C, res2);
		//cout << "db0\n";
		if (at1_h && at2_h
				&& at1_n && at2_n
				&& at1_o && at2_o
				&& at1_c && at2_c){
			//cout << "db1\n";
			if (at1_h->isActiveAndSet() && at2_h->isActiveAndSet()
					&& at1_n->isActiveAndSet() && at2_n->isActiveAndSet()
					&& at1_o->isActiveAndSet() && at2_o->isActiveAndSet()
					&& at1_c->isActiveAndSet() && at2_c->isActiveAndSet()){
				//cout << "db2\n";
				//cout << "db3\n";
				/*
				vector3d n1_h1_vec = at1_h->get_coords() - at1_n->get_coords();
				vector3d n2_h2_vec = at2_h->get_coords() - at2_n->get_coords();

				vector3d c1_o1_vec = at1_o->get_coords() - at1_c->get_coords();
				vector3d c2_o2_vec = at2_o->get_coords() - at2_c->get_coords();
				 */

				//if ((n1_h1_vec.dot(c2_o2_vec) < 0) && (c1_o1_vec.dot(n2_h2_vec) < 0)){
				if (angle(at1_n->get_coords(), at1_h->get_coords(), at2_o->get_coords()) > bb_angle_cutoff
						&& angle(at1_h->get_coords(), at2_o->get_coords(), at2_c->get_coords()) > bb_angle_cutoff
						&& angle(at2_n->get_coords(), at2_h->get_coords(), at1_o->get_coords()) > bb_angle_cutoff
						&& angle(at2_h->get_coords(), at1_o->get_coords(), at1_c->get_coords()) > bb_angle_cutoff){
					//cout << "db4\n";
					//|| (n2_h2_vec.dot(c1_o1_vec) < 0 && c2_o2_vec.dot(n1_h1_vec) < 0) ){
					return this->get_pair_energy(bb_lower_cutoff, bb_upper_cutoff, at1_h, at2_o, weight)
							+ this->get_pair_energy(bb_lower_cutoff, bb_upper_cutoff, at2_h, at1_o, weight);
					/*
					const double h1_o2_dist = (at1_h->get_coords() - at2_o->get_coords()).mod();
					const double h2_o1_dist = (at2_h->get_coords() - at1_o->get_coords()).mod();
					const double enrg = (weight * get_dist_energy(bb_lower_cutoff, bb_upper_cutoff, h1_o2_dist))
							+ (weight * get_dist_energy(bb_lower_cutoff, bb_upper_cutoff, h2_o1_dist)) ;
					return enrg;
					 */
				}
				else {
					//cout << "db5\n";
					return penalty; // + this->get_ca_pair_energy(pose_, res1, res2, is_parallel, 1.0 * weight);
				}


			} //
		}
	}
	else {
		//PARALLEL STRANDS
		//cout << "parallel\t" << weight << "\n";
		atom_shared_ptr at1_h = pose_->get_bb_atom(POSE::H, res1);
		atom_shared_ptr at2_h = pose_->get_bb_atom(POSE::H, res2);
		atom_shared_ptr at1_n = pose_->get_bb_atom(POSE::N, res1);
		atom_shared_ptr at2_n = pose_->get_bb_atom(POSE::N, res2);
		atom_shared_ptr at1_o = pose_->get_bb_atom(POSE::O, res1);
		atom_shared_ptr at2_o = pose_->get_bb_atom(POSE::O, res2);
		atom_shared_ptr at1_c = pose_->get_bb_atom(POSE::C, res1);
		atom_shared_ptr at2_c = pose_->get_bb_atom(POSE::C, res2);

		const const_residue_shared_ptr res1_ptr = pose_->get_residue(res1);
		const const_residue_shared_ptr res2_ptr = pose_->get_residue(res2);

		const const_residue_shared_ptr res1_m1_ptr = res1_ptr->get_prev_residue();
		//bool res1_m1OK = false;
		atom_shared_ptr at1_m1_o = atom_shared_ptr(), at1_m1_c = atom_shared_ptr();
		if (res1_m1_ptr){
				//res1_m1OK = true;
				at1_m1_o = pose_->get_bb_atom(POSE::O, res1-1);
				at1_m1_c = pose_->get_bb_atom(POSE::C, res1-1);
		}

		const const_residue_shared_ptr res1_p1_ptr = res1_ptr->get_next_residue();
		//bool res1_p1OK = false;
		atom_shared_ptr at1_p1_h = atom_shared_ptr(), at1_p1_n = atom_shared_ptr();
		if (res1_p1_ptr){
				//res1_p1OK = true;
				at1_p1_h = pose_->get_bb_atom(POSE::H, res1+1);
				at1_p1_n = pose_->get_bb_atom(POSE::N, res1+1);
		}

		const const_residue_shared_ptr res2_m1_ptr = res2_ptr->get_prev_residue();
		//bool res2_m1OK = false;
		atom_shared_ptr at2_m1_o = atom_shared_ptr(), at2_m1_c = atom_shared_ptr();
		if (res2_m1_ptr){
				//res2_m1OK = true;
				at2_m1_o = pose_->get_bb_atom(POSE::O, res2-1);
				at2_m1_c = pose_->get_bb_atom(POSE::C, res2-1);
		}

		const const_residue_shared_ptr res2_p1_ptr = res2_ptr->get_next_residue();
		//bool res2_p1OK = false;
		atom_shared_ptr at2_p1_h = atom_shared_ptr(), at2_p1_n = atom_shared_ptr();
		if (res2_p1_ptr){
				//res2_p1OK = true;
				at2_p1_h = pose_->get_bb_atom(POSE::H, res2+1);
				at2_p1_n = pose_->get_bb_atom(POSE::N, res2+1);
		}

		/*
		 * Check for at1_h -> at2_m1_o && at1_o -> at2_p1_h H-bonding
		 */
		if (at1_h->isActiveAndSet() && at2_m1_o->isActiveAndSet() && at1_o->isActiveAndSet() && at2_p1_h->isActiveAndSet()
				&& at2_m1_c->isActiveAndSet() && at1_n->isActiveAndSet() && at2_p1_n->isActiveAndSet() && at1_c->isActiveAndSet()
				&& angle(at1_h->get_coords(), at2_m1_o->get_coords(), at2_m1_c->get_coords()) > bb_angle_cutoff
				&& angle(at1_n->get_coords(), at1_h->get_coords(), at2_m1_o->get_coords()) > bb_angle_cutoff
				&& angle(at1_o->get_coords(), at2_p1_h->get_coords(), at2_p1_n->get_coords()) > bb_angle_cutoff
				&& angle(at1_c->get_coords(), at1_o->get_coords(), at2_p1_h->get_coords()) > bb_angle_cutoff){
			const double val = this->get_pair_energy(bb_lower_cutoff, bb_upper_cutoff, at1_h, at2_m1_o, weight)
					+ this->get_pair_energy(bb_lower_cutoff, bb_upper_cutoff, at1_o, at2_p1_h, weight);
			//cout << "db0\t" << val << "\n";
			return val;
		}
		/*
		 * Check for at2_h -> at1_m1_o && at2_o -> at1_p1_h H-bonding
		 */
		else if (at2_h->isActiveAndSet() && at1_m1_o->isActiveAndSet() && at2_o->isActiveAndSet() && at1_p1_h->isActiveAndSet()
				&& at1_m1_c->isActiveAndSet() && at2_n->isActiveAndSet() && at1_p1_n->isActiveAndSet() && at2_c->isActiveAndSet()
				&& angle(at2_h->get_coords(), at1_m1_o->get_coords(), at1_m1_c->get_coords()) > bb_angle_cutoff
				&& angle(at2_n->get_coords(), at2_h->get_coords(), at1_m1_o->get_coords()) > bb_angle_cutoff
				&& angle(at2_o->get_coords(), at1_p1_h->get_coords(), at1_p1_n->get_coords()) > bb_angle_cutoff
				&& angle(at2_c->get_coords(), at2_o->get_coords(), at1_p1_h->get_coords()) > bb_angle_cutoff){
			const double val = this->get_pair_energy(bb_lower_cutoff, bb_upper_cutoff, at2_h, at1_m1_o, weight)
					+ this->get_pair_energy(bb_lower_cutoff, bb_upper_cutoff, at2_o, at1_p1_h, weight);
			//cout << "db1\t" << val << "\n";
			return val;
		}
		else {
			const double val = penalty; // + this->get_ca_pair_energy(pose_, res1, res2, is_parallel, 1.0 * weight);
			//cout << "db2\t" << val << "\n";
			return val;
		}

	}
	return 0;
}

double strand_pairing_rst::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	const PRODART::POSE::META::int_int_tup_bool_double_tup_map& rsts = pose_meta_->get_strand_residue_pair_restraints();

	pose_shared_ptr pose_ = pose_meta_->get_pose();

	double total_energy = 0;


	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		for (int_int_tup_bool_double_tup_map::const_iterator it = rsts.begin(); it != rsts.end(); it++){
			const int res1 = it->first.get<0>();
			const int res2 = it->first.get<1>();
			const bool is_parallel = it->second.get<0>();
			const double weight = it->second.get<1>();
			total_energy += get_ca_pair_energy(pose_, res1, res2, is_parallel, weight);
		}
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta){
		for (int_int_tup_bool_double_tup_map::const_iterator it = rsts.begin(); it != rsts.end(); it++){
			//cout << "getting energy..." << endl;
			const int res1 = it->first.get<0>();
			const int res2 = it->first.get<1>();
			const bool is_parallel = it->second.get<0>();
			const double weight = it->second.get<1>();
			const double this_energy = get_bb_pair_energy(pose_, res1, res2, is_parallel, weight)
					+ get_ca_pair_energy(pose_, res1, res2, is_parallel, weight)
					+ get_cb_pair_energy(pose_, res1, res2, is_parallel, weight);
			//cout << res1 << " " << res2 << " " << this_energy << "\n";
			total_energy += this_energy;
		}
	}


	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double strand_pairing_rst::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const PRODART::POSE::META::int_int_tup_bool_double_tup_map& rsts = pose_meta_->get_strand_residue_pair_restraints();

	pose_shared_ptr pose_ = pose_meta_->get_pose();

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double pot_weight = energies_map.get_weight(this->name_vector[0]);

	double total_energy = 0;


	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		for (int_int_tup_bool_double_tup_map::const_iterator it = rsts.begin(); it != rsts.end(); it++){
			const int res1 = it->first.get<0>();
			const int res2 = it->first.get<1>();
			const bool is_parallel = it->second.get<0>();
			const double weight = it->second.get<1>();
			total_energy += get_ca_pair_energy_grad(pose_, res1, res2, is_parallel, weight, grad, pot_weight);
		}
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta){
		for (int_int_tup_bool_double_tup_map::const_iterator it = rsts.begin(); it != rsts.end(); it++){
			//cout << "getting energy..." << endl;
			const int res1 = it->first.get<0>();
			const int res2 = it->first.get<1>();
			const bool is_parallel = it->second.get<0>();
			const double weight = it->second.get<1>();
			const double this_energy = get_bb_pair_energy_grad(pose_, res1, res2, is_parallel, weight, grad, pot_weight)
					+ get_ca_pair_energy_grad(pose_, res1, res2, is_parallel, weight, grad, pot_weight)
					+ get_cb_pair_energy_grad(pose_, res1, res2, is_parallel, weight, grad, pot_weight);
			//cout << res1 << " " << res2 << " " << this_energy << "\n";
			total_energy += this_energy;
		}
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}


}
}
}

