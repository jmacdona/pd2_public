/*
 * bb_old_hbond.cpp
 *
 *  Created on: 26 Oct 2010
 *      Author: jmacdona
 */

#include "bb_old_hbond.h"

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
namespace BB{

potential_shared_ptr new_bb_old_hbond(){
	potential_shared_ptr ptr(new bb_old_hbond());
	return ptr;
}

const double bb_old_hbond::dist_bin_size(0.05);
const double bb_old_hbond::psi_bin_size(5.0 * (UTILS::PI/180.0));
const double bb_old_hbond::omega_bin_size(5.0 * (UTILS::PI/180.0));
const double bb_old_hbond::tau_bin_size(5.0 * (UTILS::PI/180.0));

const double bb_old_hbond::max_dist(2.6),
		bb_old_hbond::min_dist(1.4);
const int bb_old_hbond::max_dist_bin(static_cast<int>(ceil(max_dist / dist_bin_size))),
	bb_old_hbond::min_dist_bin(static_cast<int>(ceil(min_dist / dist_bin_size))),
	bb_old_hbond::max_psi_bin(static_cast<int>(ceil(UTILS::PI / psi_bin_size))),
	bb_old_hbond::max_omega_bin(static_cast<int>(ceil(UTILS::PI / omega_bin_size))),
	bb_old_hbond::max_tau_bin(static_cast<int>(ceil((2.0 * UTILS::PI)  / tau_bin_size)));

const double bb_old_hbond::fade_off_dist(4.0);
const double bb_old_hbond::angular_fade_off_dist(2.2);
const double bb_old_hbond::zero_dist_energy(20.0);

const double bb_old_hbond::grad_h = 1e-8;

bb_old_hbond::bb_old_hbond(){
	this->name_vector.clear();
	this->name_vector.resize((int)hb_not, potentials_name("placeholder"));
	this->name_vector[(int)hb_other] = potentials_name("bb_old_hb_other");
	this->name_vector[(int)hb_helix] = potentials_name("bb_old_hb_helix");
	this->name_vector[(int)hb_strand] = potentials_name("bb_old_hb_strand");


	/*
	cout << UTILS::PI / psi_bin_size << endl;
	cout << (int)hb_not << endl;
	cout << max_dist_bin << endl;
	cout << max_psi_bin << endl;
	cout << max_omega_bin << endl;
	cout << max_tau_bin << endl;
	*/

	H_O_dist_bin_score.resize((int)hb_not, double_vector(max_dist_bin, 0));
	psi_bin_score.resize((int)hb_not, double_vector(max_psi_bin, 0));
	omega_bin_score.resize((int)hb_not, double_vector(max_omega_bin, 0));
	tau_bin_score.resize((int)hb_not, double_vector(max_tau_bin, 0));
	hbond_energy.resize((int)hb_not, 0);

}

bool bb_old_hbond::init(){

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:bb_old_hbond");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	this->load_data(inputdb);
	inputdb.close();


	return true;
}






double bb_old_hbond::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;
	double_vector enrg_components((int)hb_not , 0);

	const nb_ele_vector& pair_list = bb_meta_dat->get_H_O_pair_list();
	nb_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		if (iter->dist <= this->fade_off_dist){
			const hb_type type = this->get_hbond_type(*iter, bb_meta_dat);
			if (type != hb_not){
				enrg_components[(int)type] += this->get_energy(*iter, type, bb_meta_dat);
			}
		}
	}

	for (int i = 0; i < (int)hb_not; i++) {
		total_energy += energies_map.add_energy_component(name_vector[i], enrg_components[i]);
	}
	return total_energy;
}

double bb_old_hbond::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	double weights[(int)hb_not];
	weights[(int)hb_other] = energies_map.get_weight(this->name_vector[(int)hb_other]);
	weights[(int)hb_helix] = energies_map.get_weight(this->name_vector[(int)hb_helix]);
	weights[(int)hb_strand] = energies_map.get_weight(this->name_vector[(int)hb_strand]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;
	double_vector enrg_components((int)hb_not , 0);

	const nb_ele_vector& pair_list = bb_meta_dat->get_H_O_pair_list();
	nb_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		if (iter->dist <= this->fade_off_dist){
			const hb_type type = this->get_hbond_type(*iter, bb_meta_dat);
			if (type != hb_not) enrg_components[(int)type] += this->get_energy_grad(*iter, type, bb_meta_dat, grad, weights[(int)type]);
		}
	}

	for (int i = 0; i < (int)hb_not; i++) {
		total_energy += energies_map.add_energy_component(name_vector[i], enrg_components[i]);
	}

	//assert(is_vector3d_vector_valid(grad));

	return total_energy;//energies_map.add_energy_component(name_vector[0], total_energy);
}


bb_old_hbond::hb_type bb_old_hbond::get_hbond_type(const PRODART::POSE::META::nb_pair_element& ele,
		const bb_pose_meta_shared_ptr bb_meta_dat) const{

	PRODART::POSE::four_state_sec_struct bin_i = bb_meta_dat->get_residue_conf_class(ele.res_num1);//this->getPhiPsiBin(phi_i,psi_i);

	PRODART::POSE::four_state_sec_struct bin_j = bb_meta_dat->get_residue_conf_class(ele.res_num2);//this->getPhiPsiBin(phi_j,psi_j);

	if (ele.seq_sep >= bb_old_hbond::min_seq_sep){
		if ((bin_i == ss4_HELIX && bin_j == ss4_STRAND) || (bin_j == ss4_HELIX && bin_i == ss4_STRAND)){
			return hb_other;
		}
		else if (ele.seq_sep == 4 && (bin_i == ss4_HELIX || bin_j == ss4_HELIX)){
			return hb_helix;//(ele.sec_struct = 1);
		}
		else if (bin_i == ss4_STRAND || bin_j == ss4_STRAND ){
			return hb_strand;
		}
		else {
			return hb_other;
		}
	}

	return hb_not;
}


double bb_old_hbond::get_dist_energy(hb_type type, const double dist) const{

	//double total_energy = 0;

	const int dist_bin = this->get_dist_bin(dist);
	const double central_dist = this->get_cdist_from_bin(dist_bin);

	if ( dist <= this->get_cdist_from_bin(min_dist_bin) ){

		const double dist_energy = H_O_dist_bin_score[(int)type][min_dist_bin - 1]  + hbond_energy[(int)type];
		const double gradient = (dist_energy - zero_dist_energy) /
								( this->get_cdist_from_bin(min_dist_bin) );

		const double energy = zero_dist_energy + (dist * gradient);

		return energy;

	}
	else if ( dist < this->get_cdist_from_bin(max_dist_bin - 1) ){

		if (dist > central_dist){
			const double central_energy = H_O_dist_bin_score[(int)type][dist_bin];//other_H_O_dist_bin_score[dist_bin1];
			const double central_energy_p1 = H_O_dist_bin_score[(int)type][dist_bin+1];//other_H_O_dist_bin_score[dist_bin1 + 1];

			const double gradient = (central_energy_p1 - central_energy) / dist_bin_size;
			const double energy = (central_energy + ((dist - central_dist) * gradient)) + hbond_energy[(int)type];

			return energy;
		}
		else {
			const double central_dist_m1 = this->get_cdist_from_bin(dist_bin - 1);
			const double central_energy = H_O_dist_bin_score[(int)type][dist_bin];//other_H_O_dist_bin_score[dist_bin1];
			const double central_energy_m1 = H_O_dist_bin_score[(int)type][dist_bin-1];//other_H_O_dist_bin_score[dist_bin1 - 1];

			const double gradient = (central_energy - central_energy_m1) / dist_bin_size;
			const double energy = (central_energy_m1 + ((dist - central_dist_m1) * gradient)) + hbond_energy[(int)type];

			return energy;
		}

	}
	else if (dist < fade_off_dist){
		const double dist_energy = this->H_O_dist_bin_score[(int)type][max_dist_bin - 1]  + hbond_energy[(int)type];
		const double gradient = dist_energy /
								( fade_off_dist - this->get_cdist_from_bin(max_dist_bin- 1) );

		const double energy = ((fade_off_dist - dist) * gradient);

		return energy;
	}

	//return total_energy;
	return 0;

}

double bb_old_hbond::get_omega_energy(hb_type type, const double omega) const{
	const int omega_bin = this->get_omega_bin(omega);
	const double central_omega = this->get_comega_from_bin(omega_bin);
	double central_energy = this->omega_bin_score[(int)type][omega_bin];

	if (omega > central_omega){
		const int omega_bin_p1 = omega_bin + 1;
		if (omega_bin_p1 < max_omega_bin ){
			const double central_energy_p1 = this->omega_bin_score[(int)type][omega_bin_p1];
			const double gradient = (central_energy_p1 - central_energy) / omega_bin_size;
			const double energy = central_energy + ((omega - central_omega) * gradient);

			return energy;
		}
		else {
			const double central_energy_m1 = omega_bin_score[(int)type][omega_bin - 1];
			const double gradient = (central_energy - central_energy_m1) / omega_bin_size;
			const double energy = central_energy + ((omega - central_omega) * gradient);

			return energy;

		}
	}
	else {
		const int omega_bin_m1 = omega_bin - 1;
		if (omega_bin_m1 >= 0){
			const double central_energy_m1 = this->omega_bin_score[(int)type][omega_bin_m1];
			const double central_omega_m1 = this->get_comega_from_bin(omega_bin - 1);
			const double gradient = (central_energy - central_energy_m1) / omega_bin_size;
			const double energy = central_energy_m1 + ((omega - central_omega_m1) * gradient);

			return energy;
		}
		else{
			const double central_energy_p1 = this->omega_bin_score[(int)type][omega_bin +1];
			const double gradient = (central_energy_p1 - central_energy) / omega_bin_size;
			const double energy = central_energy + ((omega - central_omega) * gradient);

			return energy;
		}

	}
	return 0;
}
double bb_old_hbond::get_psi_energy(hb_type type, const double psi) const{
	const int psi_bin = this->get_psi_bin(psi);
	const double central_psi = this->get_cpsi_from_bin(psi_bin);
	double central_energy = this->psi_bin_score[(int)type][psi_bin];

	if (psi > central_psi){
		const int psi_bin_p1 = psi_bin + 1;
		if (psi_bin_p1 < max_psi_bin ){
			const double central_energy_p1 = this->psi_bin_score[(int)type][psi_bin_p1];
			const double gradient = (central_energy_p1 - central_energy) / psi_bin_size;
			const double energy = central_energy + ((psi - central_psi) * gradient);

			return energy;
		}
		else {
			const double central_energy_m1 = psi_bin_score[(int)type][psi_bin - 1];
			const double gradient = (central_energy - central_energy_m1) / psi_bin_size;
			const double energy = central_energy + ((psi - central_psi) * gradient);

			return energy;

		}
	}
	else {
		const int psi_bin_m1 = psi_bin - 1;
		if (psi_bin_m1 >= 0){
			const double central_energy_m1 = this->psi_bin_score[(int)type][psi_bin_m1];
			const double central_psi_m1 = this->get_cpsi_from_bin(psi_bin - 1);
			const double gradient = (central_energy - central_energy_m1) / psi_bin_size;
			const double energy = central_energy_m1 + ((psi - central_psi_m1) * gradient);

			return energy;
		}
		else{
			const double central_energy_p1 = this->psi_bin_score[(int)type][psi_bin +1];
			const double gradient = (central_energy_p1 - central_energy) / psi_bin_size;
			const double energy = central_energy + ((psi - central_psi) * gradient);

			return energy;
		}

	}
}
double bb_old_hbond::get_tau_energy(hb_type type, const double tau) const{
	const int tau_bin = this->get_tau_bin(tau);
	const double central_tau = this->get_ctau_from_bin(tau_bin);
	double central_energy = this->tau_bin_score[(int)type][tau_bin];

	if (tau > central_tau){
		const int tau_bin_p1 = tau_bin + 1;
		if (tau_bin_p1 < max_tau_bin ){
			const double central_energy_p1 = this->tau_bin_score[(int)type][tau_bin_p1];
			const double gradient = (central_energy_p1 - central_energy) / tau_bin_size;
			const double energy = central_energy + ((tau - central_tau) * gradient);

			return energy;
		}
		else {
			const double central_energy_p1 = this->tau_bin_score[(int)type][0];
			const double gradient = (central_energy_p1 - central_energy) / tau_bin_size;
			const double energy = central_energy + ((tau - central_tau) * gradient);

			return energy;

		}
	}
	else {
		const int tau_bin_m1 = tau_bin - 1;
		if (tau_bin_m1 >= 0){
			const double central_energy_m1 = this->tau_bin_score[(int)type][tau_bin_m1];
			const double central_tau_m1 = this->get_ctau_from_bin(tau_bin - 1);
			const double gradient = (central_energy - central_energy_m1) / tau_bin_size;
			const double energy = central_energy_m1 + ((tau - central_tau_m1) * gradient);

			return energy;
		}
		else{
			const double central_energy_m1 = this->tau_bin_score[(int)type][max_tau_bin - 1];
			const double central_tau_m1 = this->get_ctau_from_bin(tau_bin - 1);
			const double gradient = (central_energy - central_energy_m1) / tau_bin_size;
			const double energy = central_energy_m1 + ((tau - central_tau_m1) * gradient);

			return energy;
		}

	}
}

double bb_old_hbond::get_angle_fade_off(const double dist) const{

	if (dist >= this->get_cdist_from_bin(max_dist_bin - 1)){
		return 0;
	}
	else if (dist > angular_fade_off_dist){
		const double gradient = -1.0
			/ ( this->get_cdist_from_bin(max_dist_bin - 1) - angular_fade_off_dist );
		return 1.0 + ((dist - angular_fade_off_dist) * gradient);
	}
	else {
		return 1.0;
	}
}

double bb_old_hbond::get_energy(const hb_type type,
		const double dist,
		const double omega,
		const double psi,
		const double tau) const{

	const double dist_energy = get_dist_energy(type, dist);
	const double omega_energy = this->get_omega_energy(type, omega);
	const double psi_energy = this->get_psi_energy(type, psi);
	const double tau_energy = this->get_tau_energy(type, tau);

	const double angle_energy = omega_energy
			+ psi_energy
			+ tau_energy;

	const double angle_weight = this->get_angle_fade_off(dist);

	/*
	cout << "debug_weight\t" << ele.dist << "\t" << angle_weight << "\t" << hb_type_to_string(type) << endl;
	cout << "debug_dist\t" << ele.dist << "\t" << dist_energy << "\t" << hb_type_to_string(type) << endl;
	cout << "debug_omega\t" << omega << "\t" << omega_energy << "\t" << hb_type_to_string(type) << endl;
	cout << "debug_psi\t" << psi << "\t" << psi_energy << "\t" << hb_type_to_string(type) << endl;
	cout << "debug_tau\t" << tau << "\t" << tau_energy << "\t" << hb_type_to_string(type) << endl;
	*/

	return (angle_weight * angle_energy) + dist_energy;
}

double bb_old_hbond::get_energy(const hb_type type,
		const double dist,
		const UTILS::vector3d n,
		const UTILS::vector3d h,
		const UTILS::vector3d ca,
		const UTILS::vector3d c,
		const UTILS::vector3d o) const{

	const double omega = angle(n, h, o);
	const double psi = angle(h, o, c);
	const double tau = dihedral(h, o, c, ca);

	return get_energy(type,
			dist,
			omega,
			psi,
			tau);

}

double bb_old_hbond::get_energy(const PRODART::POSE::META::nb_pair_element& ele,
		hb_type type,
		const PRODART::POSE::META::bb_pose_meta_shared_ptr bb_meta_dat) const{

	double total_energy = 0;

	if (ele.dist < fade_off_dist){
		if ((ele.atype1 == bb_pose_meta::bb_aty_O() && ele.atype2 == bb_pose_meta::bb_aty_H())
				|| (ele.atype2 == bb_pose_meta::bb_aty_O() && ele.atype1 == bb_pose_meta::bb_aty_H())){
			const_atom_shared_ptr n;
			const_atom_shared_ptr h;
			const_atom_shared_ptr ca;
			const_atom_shared_ptr c;
			const_atom_shared_ptr o;

			if (ele.atype1 == bb_pose_meta::bb_aty_O() && ele.atype2 == bb_pose_meta::bb_aty_H()){
				h = ele.atom2_ptr;
				o = ele.atom1_ptr;

				n = bb_meta_dat->get_pose()->get_bb_atom(POSE::N, ele.res_num2);

				ca = bb_meta_dat->get_pose()->get_bb_atom(POSE::CA, ele.res_num1);
				c = bb_meta_dat->get_pose()->get_bb_atom(POSE::C, ele.res_num1);

			}
			else {
				h = ele.atom1_ptr;
				o = ele.atom2_ptr;

				n = bb_meta_dat->get_pose()->get_bb_atom(POSE::N, ele.res_num1);

				ca = bb_meta_dat->get_pose()->get_bb_atom(POSE::CA, ele.res_num2);
				c = bb_meta_dat->get_pose()->get_bb_atom(POSE::C, ele.res_num2);

			}
			if (n && h && ca && c && o){
				if (n->isActive() && h->isActive() && ca->isActive() && c->isActive() && o->isActive()){


					/*
					const double omega = angle(n->get_coords(), h->get_coords(), o->get_coords());
					const double psi = angle(h->get_coords(), o->get_coords(), c->get_coords());
					const double tau = dihedral(h->get_coords(), o->get_coords(), c->get_coords(), ca->get_coords());
					*/

					/*
					const double dist_energy = get_dist_energy(type, ele.dist);
					const double omega_energy = this->get_omega_energy(type, omega);
					const double psi_energy = this->get_psi_energy(type, psi);
					const double tau_energy = this->get_tau_energy(type, tau);

					const double angle_energy = omega_energy
							+ psi_energy
							+ tau_energy;

					const double angle_weight = this->get_angle_fade_off(ele.dist);
					*/

					/*
					cout << "debug_weight\t" << ele.dist << "\t" << angle_weight << "\t" << hb_type_to_string(type) << endl;
					cout << "debug_dist\t" << ele.dist << "\t" << dist_energy << "\t" << hb_type_to_string(type) << endl;
					cout << "debug_omega\t" << omega << "\t" << omega_energy << "\t" << hb_type_to_string(type) << endl;
					cout << "debug_psi\t" << psi << "\t" << psi_energy << "\t" << hb_type_to_string(type) << endl;
					cout << "debug_tau\t" << tau << "\t" << tau_energy << "\t" << hb_type_to_string(type) << endl;
					*/

					total_energy += get_energy(type,
							ele.dist,
							n->get_coords(),
							h->get_coords(),
							ca->get_coords(),
							c->get_coords(),
							o->get_coords());

							/*
							get_energy(type,
							ele.dist,
							omega,
							psi,
							tau);//(angle_weight * angle_energy) + dist_energy;
							*/

				}
			}

		}
	}

	return total_energy;
}

double bb_old_hbond::get_energy_grad(const PRODART::POSE::META::nb_pair_element& ele,
		hb_type type,
		const PRODART::POSE::META::bb_pose_meta_shared_ptr bb_meta_dat,
		UTILS::vector3d_vector& grad,
		const double weight) const{


	double central_value = 0;

	if (ele.dist < fade_off_dist){
		if ((ele.atype1 == bb_pose_meta::bb_aty_O() && ele.atype2 == bb_pose_meta::bb_aty_H())
				|| (ele.atype2 == bb_pose_meta::bb_aty_O() && ele.atype1 == bb_pose_meta::bb_aty_H())){
			const_atom_shared_ptr n;
			const_atom_shared_ptr h;
			const_atom_shared_ptr ca;
			const_atom_shared_ptr c;
			const_atom_shared_ptr o;

			if (ele.atype1 == bb_pose_meta::bb_aty_O() && ele.atype2 == bb_pose_meta::bb_aty_H()){
				h = ele.atom2_ptr;
				o = ele.atom1_ptr;

				n = bb_meta_dat->get_pose()->get_bb_atom(POSE::N, ele.res_num2);

				ca = bb_meta_dat->get_pose()->get_bb_atom(POSE::CA, ele.res_num1);
				c = bb_meta_dat->get_pose()->get_bb_atom(POSE::C, ele.res_num1);

			}
			else {
				h = ele.atom1_ptr;
				o = ele.atom2_ptr;

				n = bb_meta_dat->get_pose()->get_bb_atom(POSE::N, ele.res_num1);

				ca = bb_meta_dat->get_pose()->get_bb_atom(POSE::CA, ele.res_num2);
				c = bb_meta_dat->get_pose()->get_bb_atom(POSE::C, ele.res_num2);

			}
			if (n && h && ca && c && o){
				if (n->isActive() && h->isActive() && ca->isActive() && c->isActive() && o->isActive()){

					/*
					const double central_omega = angle(n->get_coords(), h->get_coords(), o->get_coords());
					const double central_psi = angle(h->get_coords(), o->get_coords(), c->get_coords());
					const double central_tau = dihedral(h->get_coords(), o->get_coords(), c->get_coords(), ca->get_coords());
					*/

					central_value += get_energy(type,
							ele.dist,
							n->get_coords(),
							h->get_coords(),
							ca->get_coords(),
							c->get_coords(),
							o->get_coords());



					const int indices[5] = {n->get_seq_num(),
							h->get_seq_num(),
							ca->get_seq_num(),
							c->get_seq_num(),
							o->get_seq_num()};
					const vector3d init_vec[5] = {n->get_coords(),
							h->get_coords(),
							ca->get_coords(),
							c->get_coords(),
							o->get_coords()};
					vector3d vec[5] = {n->get_coords(),
							h->get_coords(),
							ca->get_coords(),
							c->get_coords(),
							o->get_coords()};

					for (int v = 0; v < 5; v++){
						for (int i = 0; i <3; i++){
							double temp = init_vec[v][i] + grad_h;
							dummy_funct(temp);
							const double hh =  temp - vec[v][i];

							vec[v][i] = init_vec[v][i];
							vec[v][i] += hh;
							double dist = (vec[1] - vec[4]).mod();
							const double new_val = get_energy(type, dist, vec[0], vec[1], vec[2], vec[3], vec[4]);//this->get_energy(ele, new_dih);

							vec[v][i] = init_vec[v][i];
							vec[v][i] -= hh;
							dist = (vec[1] - vec[4]).mod();
							const double new_low_val = get_energy(type, dist, vec[0], vec[1], vec[2], vec[3], vec[4]);//this->get_energy(ele, new_low_dih);

							vec[v][i] = init_vec[v][i];
							const double diff = new_val - new_low_val;
							const double this_grad = diff / (2.0 * hh);
							(grad[indices[v]])[i] += weight * this_grad;
						}
					}

				}
			}

		}
	}

	return central_value;

}



int bb_old_hbond::get_dist_bin(double val) const{
	return static_cast<int>(val / bb_old_hbond::dist_bin_size);
}

int bb_old_hbond::get_omega_bin(double val) const{
	return static_cast<int>(val / bb_old_hbond::omega_bin_size);
}

int bb_old_hbond::get_psi_bin(double val) const{
	return static_cast<int>(val / bb_old_hbond::psi_bin_size);
}

int bb_old_hbond::get_tau_bin(double val) const{
	return static_cast<int>( (val + PI) / bb_old_hbond::tau_bin_size);
}


double bb_old_hbond::get_cdist_from_bin(int bin) const{
	return (bb_old_hbond::dist_bin_size * static_cast<int>(bin)) + (0.5 * bb_old_hbond::dist_bin_size);
}

double bb_old_hbond::get_comega_from_bin(int bin) const{
	return (bb_old_hbond::omega_bin_size * static_cast<int>(bin)) + (0.5 * bb_old_hbond::omega_bin_size);
}
double bb_old_hbond::get_cpsi_from_bin(int bin) const{
	return (bb_old_hbond::psi_bin_size * static_cast<int>(bin)) + (0.5 * bb_old_hbond::psi_bin_size);
}
double bb_old_hbond::get_ctau_from_bin(int bin) const{
	return ((bb_old_hbond::tau_bin_size * static_cast<int>(bin)) + (0.5 * bb_old_hbond::tau_bin_size)) - PI;
}

double bb_old_hbond::get_min_dist_from_bin(int bin) const{
	return bb_old_hbond::dist_bin_size * static_cast<int>(bin);
}
double bb_old_hbond::get_min_omega_from_bin(int bin) const{
	return bb_old_hbond::omega_bin_size * static_cast<int>(bin);
}
double bb_old_hbond::get_min_psi_from_bin(int bin) const{
	return bb_old_hbond::psi_bin_size * static_cast<int>(bin);
}
double bb_old_hbond::get_min_tau_from_bin(int bin) const{
	return (bb_old_hbond::tau_bin_size * static_cast<int>(bin)) - PI;
}











std::istream& bb_old_hbond::load_data( std::istream& input ){

	string lineStr;


	long length, lineNum = 0 ;


	string_vector SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
					&& SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 7 ){

				string paraName = SplitVec[1];
				trim(paraName);

				const int bin = lexical_cast<int>(SplitVec[2]);

				//cout << bin << endl;

				if ( paraName.compare("H_O_dist") == 0 ) {
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double beta = lexical_cast<double>(SplitVec[6]);

					this->H_O_dist_bin_score[(int)hb_other][bin] = other;
					this->H_O_dist_bin_score[(int)hb_helix][bin] = alpha;
					this->H_O_dist_bin_score[(int)hb_strand][bin] = beta;

				}
				else if ( paraName.compare("psi") == 0 ){
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double beta = lexical_cast<double>(SplitVec[6]);


					this->psi_bin_score[(int)hb_other][bin] = other;
					this->psi_bin_score[(int)hb_helix][bin] = alpha;
					this->psi_bin_score[(int)hb_strand][bin] = beta;


				}
				else if ( paraName.compare("omega") == 0 ){
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double beta = lexical_cast<double>(SplitVec[6]);

					this->omega_bin_score[(int)hb_other][bin] = other;
					this->omega_bin_score[(int)hb_helix][bin] = alpha;
					this->omega_bin_score[(int)hb_strand][bin] = beta;

				}
				else if ( paraName.compare("tau") == 0 ){
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double beta = lexical_cast<double>(SplitVec[6]);

					this->tau_bin_score[(int)hb_other][bin] = other;
					this->tau_bin_score[(int)hb_helix][bin] = alpha;
					this->tau_bin_score[(int)hb_strand][bin] = beta;

				}
				else if (paraName.compare("hbond_energies") == 0){
					const double other = lexical_cast<double>(SplitVec[4]);
					const double alpha = lexical_cast<double>(SplitVec[5]);
					const double beta = lexical_cast<double>(SplitVec[6]);

					this->hbond_energy[(int)hb_other] = other;
					this->hbond_energy[(int)hb_helix] = alpha;
					this->hbond_energy[(int)hb_strand] = beta;


				}
				else {
					cout << "ERROR - unknown parameter name: " << paraName << endl;
					//fatalError = true;
					//cout << paraValue << endl;
				}


			}
		}
	}


	/*
	for (int i = 0; i < H_O_dist_bin_score.size() ; i++){
		for (int j = 0; j < H_O_dist_bin_score[i].size(); j++){
			cout << i << "\t" << j << "\t" << H_O_dist_bin_score[i][j] << endl;
		}
	}
	*/

	return input;
}


void bb_old_hbond::dummy_funct(double val) const{
	//double temp = val;
}



}
}
}
}


