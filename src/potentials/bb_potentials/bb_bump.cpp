/*
 * bb_bump.cpp
 *
 *  Created on: 26 Oct 2010
 *      Author: jmacdona
 */



#include "bb_bump.h"

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

const double tol = 1e-5;//0.0001;

potential_shared_ptr new_bb_bump(){
	potential_shared_ptr ptr(new bb_bump());
	return ptr;
}

void bb_bump::add_min_dist(const POSE::atom_type at1, const POSE::atom_type at2, const double dist ){
	atom_type2_tuple tup = boost::make_tuple(at1, at2);
	min_dists_map[tup] = dist;
	min_dists_sq_map[tup] = dist*dist;

	string key(at1.get_label());
	key.append(at2.get_label());
	min_dists_un_map[key] = dist;
	min_dists_sq_un_map[key] = dist*dist;
}

bb_bump::bb_bump() : min_bond_sep(6), vdw_deflator(0.99), default_min_dist(1.5), default_min_dist_sq(default_min_dist*default_min_dist), outer_bound(4.0) {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("bb_bump"));

	atomVdwParams[atom_type("C")] = vdw_deflator * 3.750;
	atomVdwParams[atom_type("O")] = vdw_deflator * 2.960;
	atomVdwParams[atom_type("N")] = vdw_deflator * 3.250;
	atomVdwParams[atom_type("H")] = vdw_deflator * 0.000;
	atomVdwParams[atom_type("CA")] = 3.800;
	atomVdwParams[atom_type("CB")] = vdw_deflator * 3.910;

	//min_dists_un_map[""];

	//atom_type_double_map::const_iterator it1, it2;
	atom_type_double_unord_map::const_iterator it1, it2;
	for (it1 = atomVdwParams.begin(); it1 != atomVdwParams.end(); it1++){
		for (it2 = atomVdwParams.begin(); it2 != atomVdwParams.end(); it2++){
			//atom_type2_tuple tup = boost::make_tuple(it1->first, it2->first);
			const double min_dist = sqrt(it1->second * it2->second);
			//cout << min_dist << endl;
			this->add_min_dist(it1->first, it2->first, min_dist);
			//min_dists_map[tup] = min_dist;
			//min_dists_sq_map[tup] = min_dist*min_dist;
		}
	}

	// due to clashing in h-bonding between strands
	add_min_dist(atom_type("N"), atom_type("O"), 2.81);
	add_min_dist(atom_type("O"), atom_type("N"), 2.81);

	//due to clashing in h-bonding between strands
	add_min_dist(atom_type("CA"), atom_type("O"), 3.11);
	add_min_dist(atom_type("O"), atom_type("CA"), 3.11);

	//due to clashing in local 1-5 structure in strands and in alpha helices
	add_min_dist(atom_type("C"), atom_type("O"), 2.75);
	add_min_dist(atom_type("O"), atom_type("C"), 2.75);

	//due to clashing in local 1-5 structure in strands especially PRO
	add_min_dist(atom_type("CA"), atom_type("C"), 3.725);
	add_min_dist(atom_type("C"), atom_type("CA"), 3.725);

	//due to clashing in distant structure
	add_min_dist(atom_type("CB"), atom_type("O"), 3.170);
	add_min_dist(atom_type("O"), atom_type("CB"), 3.170);

	//due to clashing in local (seq_sep 2) structure
	add_min_dist(atom_type("N"), atom_type("C"), 3.15);
	add_min_dist(atom_type("C"), atom_type("N"), 3.15);


}

bool bb_bump::init(){

	return true;
}

double bb_bump::get_min_dist(const POSE::atom_type at1, const POSE::atom_type at2) const{
	atom_type2_tuple tup = boost::make_tuple(at1, at2);
	return min_dists_map.find(tup) != min_dists_map.end() ? min_dists_map.find(tup)->second : default_min_dist;

	//string key(at1.get_label());
	//key.append(at2.get_label());
	//return min_dists_un_map.find(key) != min_dists_un_map.end() ? min_dists_un_map.find(key)->second : default_min_dist;
}
double bb_bump::get_min_dist_sq(const POSE::atom_type at1, const POSE::atom_type at2) const{
	atom_type2_tuple tup = boost::make_tuple(at1, at2);
	return min_dists_sq_map.find(tup) != min_dists_sq_map.end() ? min_dists_sq_map.find(tup)->second : default_min_dist_sq;

	//string key(at1.get_label());
	//key.append(at2.get_label());
	//return min_dists_sq_un_map.find(key) != min_dists_sq_un_map.end() ? min_dists_sq_un_map.find(key)->second : default_min_dist;
}

double bb_bump::get_energy(const PRODART::POSE::META::nb_pair_element& ele) const{

	if (ele.dist < this->outer_bound){
		const double dist_sq = ele.dist_sq;
		const double radius = get_min_dist(ele.atype1, ele.atype2);//3;
		const double radius_sq = radius*radius;//3;

		const int seq_sep = ele.seq_sep;	//99999;
		/*
	cout << ele.atype1.get_label() << "\t" << ele.atype2.get_label() << "\t" << seq_sep << "\t"
			<< ele.dist << "\t"
			<< radius << endl;
		 */
		if (dist_sq < radius_sq
				&& seq_sep >= this->min_bond_sep){
			/*
		cout << ele.atype1.get_label() << "\t" << ele.atype2.get_label() << "\t" << seq_sep << "\t"
				<< ele.dist << "\t"
				<< radius << endl;
			 */
			return (radius_sq - (dist_sq)) / radius;
			//total_score += pow((CA_radius_sq - (dist_sq)),2) / CA_radius;
		}
	}
	return 0;
}

double bb_bump::get_energy_grad(const PRODART::POSE::META::nb_pair_element& ele, UTILS::vector3d_vector& grad, const double weight ) const{
	if (ele.dist < this->outer_bound){
		const double dist_sq = ele.dist_sq;
		const int seq_sep = ele.seq_sep;	//99999;
		const double radius = get_min_dist(ele.atype1, ele.atype2);//3;
		const double radius_sq = radius*radius;//3;

		if (dist_sq < radius_sq
				&& seq_sep >= this->min_bond_sep){
			const double energy = (radius_sq - (dist_sq)) / radius;
			const double grad_mag = weight * ((2.0 * ele.dist) / radius);

#ifndef NDEBUG
			if (ele.dist <= tol){
				cerr << "ERROR: distance betweem atoms is less than tol: " << ele.dist << endl;
				cerr << ele.atom1_ptr->get_residue()->get_pdb_residue_index() << "\t"
						<< ele.atom1_ptr->get_chain()->getChainID() << "\t"
						<< ele.atom1_ptr->get_type().get_label() << "\t"
						<< ele.atom2_ptr->get_residue()->get_pdb_residue_index() << "\t"
						<< ele.atom2_ptr->get_chain()->getChainID() << "\t"
						<< ele.atom2_ptr->get_type().get_label() << "\t"
						<< endl;
			}

#endif




			const vector3d unit_grad_vec = ele.dist > tol ? (ele.atom1_ptr->get_coords() - ele.atom2_ptr->get_coords()) / ele.dist : vector3d(0,0,0);
			const vector3d atm1_grad = -grad_mag * unit_grad_vec;
			const vector3d atm2_grad = grad_mag * unit_grad_vec;

			const int index1 = ele.atom1_ptr->get_seq_num();
			const int index2 = ele.atom2_ptr->get_seq_num();

			grad[index1] += atm1_grad;
			grad[index2] += atm2_grad;

			return energy;
		}
		else{
			return 0;
		}
	}
	return 0;
}



double bb_bump::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	const const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;


	const nb_ele_vector& pair_list = bb_meta_dat->get_all_bb_pair_list();
	nb_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy(*iter);
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double bb_bump::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const bb_pose_meta_shared_ptr bb_meta_dat = static_pointer_cast<bb_pose_meta, pose_meta_interface>(pose_meta_);
	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double weight = energies_map.get_weight(this->name_vector[0]);
	const PRODART::POSE::const_pose_shared_ptr pose_ = pose_meta_->get_pose();
	double  total_energy = 0;

	const nb_ele_vector& pair_list = bb_meta_dat->get_all_bb_pair_list();
	nb_ele_vector::const_iterator iter;
	for (iter = pair_list.begin(); iter != pair_list.end(); iter++){
		total_energy += this->get_energy_grad(*iter, grad, weight);
	}

	//assert(is_vector3d_vector_valid(grad));

	return energies_map.add_energy_component(name_vector[0], total_energy);
}








}
}
}
}

