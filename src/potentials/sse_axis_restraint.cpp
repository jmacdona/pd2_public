/*
 * sse_axis_restraint.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: jmacdona
 */

#include "sse_axis_restraint.h"

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



potential_shared_ptr new_sse_axis_restraint(){
	potential_shared_ptr ptr(new sse_axis_restraint());
	return ptr;
}

const string H_pot_name("sse_helix_axis_rst");
const string E_pot_name("sse_strand_axis_rst");

sse_axis_restraint::sse_axis_restraint() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name(H_pot_name));
	this->name_vector.push_back(potentials_name(E_pot_name));
}

const double sse_axis_restraint::h  = 1e-8;

bool sse_axis_restraint::init(){

	return true;
}

double sse_axis_restraint::get_energy(const META::sse_axis_element & ele, UTILS::vector3d_vector ca_coords) const{
	POSE::four_state_sec_struct secs = ele.secs;
	vector3d p_start, p_end;
	if (secs == ss4_HELIX){
		if (!get_helix_end_points_long(ca_coords, p_start, p_end)) return 0;
	}
	else if (secs == ss4_STRAND) {
		if (!get_strand_end_points_long(ca_coords, p_start, p_end)) return 0;
	}
	else {
		return 0;
	}
	vector3d pa;
	vector3d pb;
	double mua;
	double mub;
	double  dih_ang;
	double p1_pa_pb_ang;
	double  pa_pb_p3_ang;
	double dist_sq;

	double len = (double)ca_coords.size();

	if (segmentIntersect(ele.start.coord, ele.end.coord, p_start, p_end,
			pa, pb,
			mua, mub,
			dih_ang,
			p1_pa_pb_ang,
			pa_pb_p3_ang,
			dist_sq)){

		//const double dist_sq = (pa - pb).mod_sq();

		/*
		cout << "ele.start.coord\t" << ele.start.coord << endl;
		cout << "ele.end.coord\t" << ele.end.coord << endl;
		cout << "p_start\t" << p_start << endl;
		cout << "p_end\t" << p_end << endl;
		cout << "dist_sq: " << dist_sq << endl;
		cout << "dih_ang: " << radians_to_degrees(dih_ang) << endl;
		cout << "p1_pa_pb_ang: " << radians_to_degrees(p1_pa_pb_ang) << endl;
		cout << "pa_pb_p3_ang: " << radians_to_degrees(pa_pb_p3_ang) << endl;

		PRINT_EXPR(dist_sq);
		PRINT_EXPR(pow(get_dihedral_distance(dih_ang,0),2));
		*/

		double energy = len * dist_sq * ele.seg_intersect_dist_half_bond_const;
		energy += len * pow(get_dihedral_distance(dih_ang,0),2) * ele.seg_intersect_dih_half_bond_const;
		energy += len * pow(p1_pa_pb_ang-degrees_to_radians(90.0),2) * ele.seg_intersect_ang_half_bond_const;
		energy += len * pow(pa_pb_p3_ang-degrees_to_radians(90.0),2) * ele.seg_intersect_ang_half_bond_const;

		return energy;

	}
	else {
		return 0;
	}

}

double sse_axis_restraint::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	//const PRODART::POSE::four_state_sec_struct_vector& conf_class = pose_meta_->get_conf_class();
	//const PRODART::POSE::three_state_sec_struct_vector& sec_struct = pose_meta_->get_sec_struct();
	pose_shared_ptr pose_ = pose_meta_->get_pose();
	const PRODART::POSE::META::sse_axis_element_vector& rsts = pose_meta_->get_sse_axis_restraint_list();

	double total_helix_energy = 0;
	double total_strand_energy = 0;

	for (unsigned int i = 0; i < rsts.size(); i++){
		const sse_axis_element & ele = rsts[i];

		POSE::four_state_sec_struct secs = ele.secs;
		if (secs == ss4_HELIX || secs == ss4_STRAND){
			bool all_valid = true;
			const int start_resnum = ele.start.res_num;
			const int end_resnum = ele.end.res_num;
			if (end_resnum - start_resnum >=3){
				//const vector3d start_vec = ele.start.coord;
				//const vector3d ent_vec = ele.end.coord;
				vector3d_vector ca_coords;
				ca_coords.reserve(end_resnum - start_resnum + 1);
				for (int r = start_resnum; r <= end_resnum; r++){
					atom_shared_ptr at = pose_->get_bb_atom(POSE::CA, r);
					if (at){
						if (at->isActiveAndSet()){
							ca_coords.push_back(at->get_coords());
						}
						else{
							all_valid = false;
						}
					}
					else {
						all_valid = false;
					}
				}
				if (all_valid){
					const double energy = this->get_energy(ele, ca_coords);
					if (secs == ss4_HELIX ){
						total_helix_energy += energy;//this->get_energy(ele, ca_coords);
					}
					else {
						total_strand_energy  += energy;//this->get_energy(ele, ca_coords);
					}
				}
			}
		}

	}

	return energies_map.add_energy_component(H_pot_name, total_helix_energy) + energies_map.add_energy_component(E_pot_name, total_strand_energy) ;
}

double sse_axis_restraint::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	const double helix_weight = energies_map.get_weight(H_pot_name);
	const double strand_weight = energies_map.get_weight(E_pot_name);

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	pose_shared_ptr pose_ = pose_meta_->get_pose();
	const PRODART::POSE::META::sse_axis_element_vector& rsts = pose_meta_->get_sse_axis_restraint_list();

	double total_helix_energy = 0;
	double total_strand_energy = 0;

	for (unsigned int i = 0; i < rsts.size(); i++){
		const sse_axis_element & ele = rsts[i];

		POSE::four_state_sec_struct secs = ele.secs;
		if (secs == ss4_HELIX || secs == ss4_STRAND){
			bool all_valid = true;
			const int start_resnum = ele.start.res_num;
			const int end_resnum = ele.end.res_num;
			if (end_resnum - start_resnum >=3){
				//const vector3d start_vec = ele.start.coord;
				//const vector3d ent_vec = ele.end.coord;
				vector3d_vector ca_coords;
				ca_coords.reserve(end_resnum - start_resnum + 1);
				atom_shared_ptr_vector atom_vec;
				atom_vec.reserve(end_resnum - start_resnum + 1);
				for (int r = start_resnum; r <= end_resnum; r++){
					atom_shared_ptr at = pose_->get_bb_atom(POSE::CA, r);
					if (at){
						if (at->isActiveAndSet()){
							ca_coords.push_back(at->get_coords());
							atom_vec.push_back(at);
						}
						else{
							all_valid = false;
						}
					}
					else {
						all_valid = false;
					}
				}
				if (all_valid){
					const double central_value = this->get_energy(ele, ca_coords);
					if (secs == ss4_HELIX ){
						total_helix_energy += central_value;//this->get_energy(ele, ca_coords);
					}
					else {
						total_strand_energy  += central_value;//this->get_energy(ele, ca_coords);
					}
					const double weight = secs == ss4_HELIX ? helix_weight : strand_weight;
					//*************************





					const vector3d_vector init_ca_coords = ca_coords;

					for (unsigned int an = 0; an < ca_coords.size(); an++){
						for (int r = 0; r <3; r++){

							//PRINT_EXPR(an);
							//PRINT_EXPR(r);

							double temp = init_ca_coords[an][r] + h;//init_vec[v][i] + h;
							dummy_funct(temp);
							const double hh =  temp - init_ca_coords[an][r];//vec[v][i];

							ca_coords[an][r] = init_ca_coords[an][r];
							ca_coords[an][r] += hh;
							const double new_val = this->get_energy(ele, ca_coords);//this->get_energy(new_angle, ele.equilibrium_angle, ele.half_angle_const);

							ca_coords[an][r] = init_ca_coords[an][r];
							ca_coords[an][r] -= hh;
							const double new_low_val = this->get_energy(ele, ca_coords);//this->get_energy(new_low_angle, ele.equilibrium_angle, ele.half_angle_const);

							ca_coords[an][r] = init_ca_coords[an][r];
							const double diff = new_val - new_low_val;
							const double this_grad = diff / (2.0 * hh);
							(grad[atom_vec[an]->get_seq_num()])[r] += weight * this_grad;

							//PRINT_EXPR(atom_vec[an]->get_seq_num());
							//PRINT_EXPR(r);
							//PRINT_EXPR(weight * this_grad);
							//cout << weight * this_grad << "\t";
						}
					}





					//*************************

				}
			}
		}

	}

	const double total_energy = energies_map.add_energy_component(H_pot_name, total_helix_energy) + energies_map.add_energy_component(E_pot_name, total_strand_energy);
	return total_energy;

}




void sse_axis_restraint::dummy_funct(double val) const{
	//double temp = val;
}



}
}
}
