/*
 * sec_struct_dih_ang_rst.cpp
 *
 *  Created on: 8 Nov 2011
 *      Author: jmacdona
 */
#include "sec_struct_dih_ang_rst.h"

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



potential_shared_ptr new_sec_struct_dih_ang_rst(){
	potential_shared_ptr ptr(new sec_struct_dih_ang_rst());
	return ptr;
}

const double sec_struct_dih_ang_rst::beta_min_phi = degrees_to_radians(-180),
		sec_struct_dih_ang_rst::beta_max_phi = degrees_to_radians(-30);
const double sec_struct_dih_ang_rst::beta_min_psi = degrees_to_radians(50),
		sec_struct_dih_ang_rst::beta_max_psi = degrees_to_radians(-170);

const double sec_struct_dih_ang_rst::alpha_min_phi = degrees_to_radians(-90),
		sec_struct_dih_ang_rst::alpha_max_phi = degrees_to_radians(-45);
const double sec_struct_dih_ang_rst::alpha_min_psi = degrees_to_radians(-60),
		sec_struct_dih_ang_rst::alpha_max_psi = degrees_to_radians(0);

// ideal = -151.00
const double sec_struct_dih_ang_rst::beta_min_ca_dih = degrees_to_radians(180),
		sec_struct_dih_ang_rst::beta_max_ca_dih = degrees_to_radians(-130);

// ideal = 49.7
const double sec_struct_dih_ang_rst::alpha_min_ca_dih = degrees_to_radians(48.0), //45.0
		sec_struct_dih_ang_rst::alpha_max_ca_dih = degrees_to_radians(52.0);      //55.0

sec_struct_dih_ang_rst::sec_struct_dih_ang_rst() {
	this->name_vector.clear();
	this->name_vector.push_back(potentials_name("secs_dih_ang_rst"));
}



bool sec_struct_dih_ang_rst::init(){

	return true;
}



double sec_struct_dih_ang_rst::get_dih_energy(const double equil_dih_angle, const double dih) const{
	const double energy = std::pow(UTILS::get_dihedral_distance(equil_dih_angle, dih) , 2); //ele.weight *
	return energy;
}

double sec_struct_dih_ang_rst::get_dE_ddih(const double equil_dih_angle, const double dih) const{
	const double dE = 2.0 * (UTILS::get_dihedral_distance(dih, equil_dih_angle)); // ele.weight *
	return dE;
}

double sec_struct_dih_ang_rst::get_dih_energy(POSE::atom_shared_ptr atom1_ptr,
		POSE::atom_shared_ptr atom2_ptr,
		POSE::atom_shared_ptr atom3_ptr,
		POSE::atom_shared_ptr atom4_ptr,
		const double dih,
		const double equil_ang) const{
	if (atom1_ptr
			&& atom2_ptr
			&& atom3_ptr
			&& atom4_ptr){
		if (atom1_ptr->isActive()
				&& atom2_ptr->isActive()
				&& atom3_ptr->isActive()
				&& atom4_ptr->isActive()){
			//const double dih = dihedral(ele.atom1_ptr->get_coords(), ele.atom2_ptr->get_coords(), ele.atom3_ptr->get_coords(), ele.atom4_ptr->get_coords());
			return this->get_dih_energy(equil_ang, dih);
		}
	}

		return 0;

}



double sec_struct_dih_ang_rst::get_dih_energy_ana_grad(POSE::atom_shared_ptr atom1_ptr,
		POSE::atom_shared_ptr atom2_ptr,
		POSE::atom_shared_ptr atom3_ptr,
		POSE::atom_shared_ptr atom4_ptr,
		double dih,
		const double equil_ang,
		UTILS::vector3d_vector& grad, const double weight ) const{
	if (atom1_ptr
			&& atom2_ptr
			&& atom3_ptr
			&& atom4_ptr){
		if (atom1_ptr->isActive()
				&& atom2_ptr->isActive()
				&& atom3_ptr->isActive()
				&& atom4_ptr->isActive()){

			const int indices[4] = {atom1_ptr->get_seq_num(), atom2_ptr->get_seq_num(), atom3_ptr->get_seq_num(), atom4_ptr->get_seq_num()};
			const vector3d init_vec[4] = {atom1_ptr->get_coords(), atom2_ptr->get_coords(), atom3_ptr->get_coords(), atom4_ptr->get_coords()};
			vector3d at_grad_vecs[4], f1;
			//double dih = ele.dih_angle;


			dihedral_p1_cosine_deriv(
					init_vec[0],
					init_vec[1],
					init_vec[2],
					init_vec[3],
					dih,
					//f1,
					at_grad_vecs[0]);
			dihedral_p1_cosine_deriv(
					init_vec[3],
					init_vec[2],
					init_vec[1],
					init_vec[0],
					dih,
					//f1,
					at_grad_vecs[3]);
			dihedral_p2_cosine_deriv(
					init_vec[0],
					init_vec[1],
					init_vec[2],
					init_vec[3],
					dih,
					//f1,
					at_grad_vecs[1]);
			dihedral_p2_cosine_deriv(
					init_vec[3],
					init_vec[2],
					init_vec[1],
					init_vec[0],
					dih,
					//f1,
					at_grad_vecs[2]);

			const double central_value = get_dih_energy(atom1_ptr, atom2_ptr, atom3_ptr, atom4_ptr, dih, equil_ang);
			const double dE_ddih = weight * this->get_dE_ddih(equil_ang, dih);

			/*
			cout << radians_to_degrees(ele.dih_angle) << "\t"
					<< radians_to_degrees(dihedral(init_vec[0],init_vec[1],init_vec[2],init_vec[3])) << "\t"
					<< radians_to_degrees(dih) << "\t:\t"
					<< radians_to_degrees(ele.equil_dih_angle) << "\t"
					<< radians_to_degrees(this->get_dihedral_distance(dih, ele.equil_dih_angle)) << endl;
					*/

			for (int i = 0; i < 4; i++){
				grad[indices[i]] += dE_ddih * at_grad_vecs[i];
				//cout << "ana: " << i << "\t" << dE_ddih * at_grad_vecs[i] << endl;
			}


			return central_value;
		}
	}

	return 0;
}


double sec_struct_dih_ang_rst::get_min_max_dih_energy(POSE::atom_shared_ptr atom1_ptr,
		POSE::atom_shared_ptr atom2_ptr,
		POSE::atom_shared_ptr atom3_ptr,
		POSE::atom_shared_ptr atom4_ptr,
		const double dih,
		const double min_dih,
		const double max_dih) const{
	if (is_dihedral_in_range(min_dih, max_dih, dih)){
		return 0;
	}
	else if (fabs(get_dihedral_distance(min_dih, dih)) < fabs(get_dihedral_distance(min_dih, dih))) {
		return get_dih_energy(atom1_ptr,
				atom2_ptr,
				atom3_ptr,
				atom4_ptr,
				dih,
				min_dih);
	}
	else {
		return get_dih_energy(atom1_ptr,
				atom2_ptr,
				atom3_ptr,
				atom4_ptr,
				dih,
				max_dih);
	}
}
double sec_struct_dih_ang_rst::get_min_max_dih_energy_ana_grad(POSE::atom_shared_ptr atom1_ptr,
		POSE::atom_shared_ptr atom2_ptr,
		POSE::atom_shared_ptr atom3_ptr,
		POSE::atom_shared_ptr atom4_ptr,
		double dih,
		const double min_dih,
		const double max_dih,
		UTILS::vector3d_vector& grad, const double weight ) const{
	if (is_dihedral_in_range(min_dih, max_dih, dih)){
		return 0;
	}
	else if (fabs(get_dihedral_distance(min_dih, dih)) < fabs(get_dihedral_distance(min_dih, dih))) {
		return get_dih_energy_ana_grad(atom1_ptr,
				atom2_ptr,
				atom3_ptr,
				atom4_ptr,
				dih,
				min_dih,
				grad,
				weight);
	}
	else {
		return get_dih_energy_ana_grad(atom1_ptr,
				atom2_ptr,
				atom3_ptr,
				atom4_ptr,
				dih,
				max_dih,
				grad,
				weight);
	}
}



double sec_struct_dih_ang_rst::get_energy(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{
	//const PRODART::POSE::four_state_sec_struct_vector& conf_class = pose_meta_->get_conf_class();
	//const PRODART::POSE::three_state_sec_struct_vector& sec_struct = pose_meta_->get_sec_struct();
	const PRODART::POSE::four_state_sec_struct_vector& rsts = pose_meta_->get_sec_struct_restraint_list();
	const double_vector& weights = pose_meta_->get_sec_struct_restraint_weight_list();

	double total_energy = 0;

	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		const unsigned int rescount = rsts.size();
		four_state_sec_struct prev_ss_rst = ss4_UNDEF;
		//four_state_sec_struct prev_prev_ss_rst = ss4_UNDEF;
		for (unsigned int i = 0; i < rescount; i++){
			const four_state_sec_struct ss_rst = rsts[i];
			if ( (ss_rst == ss4_STRAND || ss_rst == ss4_HELIX)
					&& ss_rst == prev_ss_rst ){
				const double weight = weights[i];
				if (i >= 1 && i+1 < rescount){
					POSE::atom_shared_ptr ca_m1 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i-1);
					POSE::atom_shared_ptr ca_0 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i);
					POSE::atom_shared_ptr ca_p1 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i+1);
					if (i+2 < rescount){
						POSE::atom_shared_ptr ca_p2 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i+2);
						if (ca_m1->get_chain() == ca_p2->get_chain() && ss_rst == rsts[i+1] && ss_rst == rsts[i+2]){
							if (ss_rst == ss4_STRAND){
								const double dih = dihedral(ca_m1->get_coords(), ca_0->get_coords(), ca_p1->get_coords(), ca_p2->get_coords());
								total_energy += weight * this->get_min_max_dih_energy(ca_m1, ca_0, ca_p1, ca_p2,
															dih, beta_min_ca_dih, beta_max_ca_dih);
							}
							else if (ss_rst == ss4_HELIX){
								const double dih = dihedral(ca_m1->get_coords(), ca_0->get_coords(), ca_p1->get_coords(), ca_p2->get_coords());
								total_energy += weight * this->get_min_max_dih_energy(ca_m1, ca_0, ca_p1, ca_p2,
															dih, alpha_min_ca_dih, alpha_max_ca_dih);
							}
						}
					}
				}
			}
			//prev_prev_ss_rst = prev_ss_rst;
			prev_ss_rst = ss_rst;
		}
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta){
		for (unsigned int i = 0; i < rsts.size(); i++){
			const four_state_sec_struct ss_rst = rsts[i];
			if (ss_rst == ss4_STRAND || ss_rst == ss4_HELIX){
				const double weight = weights[i];


				const double phi = pose_meta_->get_phi(i);
				POSE::atom_shared_ptr4_tuple phi_atoms = pose_meta_->get_pose()->get_phi_atoms(i);
				const double psi = pose_meta_->get_psi(i);
				POSE::atom_shared_ptr4_tuple psi_atoms = pose_meta_->get_pose()->get_psi_atoms(i);
				if (ss_rst == ss4_STRAND){
					total_energy += weight * this->get_min_max_dih_energy(phi_atoms.get<0>(), phi_atoms.get<1>(), phi_atoms.get<2>(), phi_atoms.get<3>(),
							phi, beta_min_phi, beta_max_phi);
					total_energy += weight * this->get_min_max_dih_energy(psi_atoms.get<0>(), psi_atoms.get<1>(), psi_atoms.get<2>(), psi_atoms.get<3>(),
							psi, beta_min_psi, beta_max_psi);
				}
				else if (ss_rst == ss4_HELIX){
					total_energy += weight * this->get_min_max_dih_energy(phi_atoms.get<0>(), phi_atoms.get<1>(), phi_atoms.get<2>(), phi_atoms.get<3>(),
							phi, alpha_min_phi, alpha_max_phi);
					total_energy += weight * this->get_min_max_dih_energy(psi_atoms.get<0>(), psi_atoms.get<1>(), psi_atoms.get<2>(), psi_atoms.get<3>(),
							psi, alpha_min_psi, alpha_max_psi);
				}

			}

		}
	}

	return energies_map.add_energy_component(name_vector[0], total_energy);
}

double sec_struct_dih_ang_rst::get_energy_with_gradient(const PRODART::POSE::META::pose_meta_shared_ptr pose_meta_,
		potentials_energies_map& energies_map) const{

	PRODART::UTILS::vector3d_vector& grad = pose_meta_->get_gradient();
	const double pot_weight = energies_map.get_weight(this->name_vector[0]);

	const PRODART::POSE::four_state_sec_struct_vector& rsts = pose_meta_->get_sec_struct_restraint_list();
	const double_vector& weights = pose_meta_->get_sec_struct_restraint_weight_list();

	double total_energy = 0;

	if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_ca_pose_meta){
		const unsigned int rescount = rsts.size();
		four_state_sec_struct prev_ss_rst = ss4_UNDEF;
		//four_state_sec_struct prev_prev_ss_rst = ss4_UNDEF;
		for (unsigned int i = 0; i < rescount; i++){
			const four_state_sec_struct ss_rst = rsts[i];
			if ( (ss_rst == ss4_STRAND || ss_rst == ss4_HELIX)
					&& ss_rst == prev_ss_rst ){
				const double weight = weights[i];
				if (i >= 1 && i+1 < rescount){
					POSE::atom_shared_ptr ca_m1 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i-1);
					POSE::atom_shared_ptr ca_0 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i);
					POSE::atom_shared_ptr ca_p1 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i+1);
					if (i+2 < rescount){
						POSE::atom_shared_ptr ca_p2 = pose_meta_->get_pose()->get_bb_atom(POSE::CA, i+2);
						if (ca_m1->get_chain() == ca_p2->get_chain() && ss_rst == rsts[i+1] && ss_rst == rsts[i+2]){
							if (ss_rst == ss4_STRAND){
								const double dih = dihedral(ca_m1->get_coords(), ca_0->get_coords(), ca_p1->get_coords(), ca_p2->get_coords());
								total_energy += weight * this->get_min_max_dih_energy_ana_grad(ca_m1, ca_0, ca_p1, ca_p2,
															dih, beta_min_ca_dih, beta_max_ca_dih,
															grad, weight * pot_weight);
							}
							else if (ss_rst == ss4_HELIX){
								const double dih = dihedral(ca_m1->get_coords(), ca_0->get_coords(), ca_p1->get_coords(), ca_p2->get_coords());
								total_energy += weight * this->get_min_max_dih_energy_ana_grad(ca_m1, ca_0, ca_p1, ca_p2,
															dih, alpha_min_ca_dih, alpha_max_ca_dih,
															grad, weight * pot_weight);
							}
						}
					}
				}
			}
			//prev_prev_ss_rst = prev_ss_rst;
			prev_ss_rst = ss_rst;
		}
	}
	else if (pose_meta_->get_pose_meta_type() == pose_meta_interface::pm_bb_pose_meta) {
		for (unsigned int i = 0; i < rsts.size(); i++){
			const four_state_sec_struct ss_rst = rsts[i];
			if (ss_rst == ss4_STRAND || ss_rst == ss4_HELIX){
				const double weight = weights[i];
				const double phi = pose_meta_->get_phi(i);
				POSE::atom_shared_ptr4_tuple phi_atoms = pose_meta_->get_pose()->get_phi_atoms(i);
				const double psi = pose_meta_->get_psi(i);
				POSE::atom_shared_ptr4_tuple psi_atoms = pose_meta_->get_pose()->get_psi_atoms(i);
				if (ss_rst == ss4_STRAND){
					total_energy += weight * this->get_min_max_dih_energy_ana_grad(phi_atoms.get<0>(), phi_atoms.get<1>(), phi_atoms.get<2>(), phi_atoms.get<3>(),
							phi, beta_min_phi, beta_max_phi,
							grad, weight * pot_weight);
					total_energy += weight * this->get_min_max_dih_energy_ana_grad(psi_atoms.get<0>(), psi_atoms.get<1>(), psi_atoms.get<2>(), psi_atoms.get<3>(),
							psi, beta_min_psi, beta_max_psi,
							grad, weight * pot_weight);
				}
				else if (ss_rst == ss4_HELIX){
					total_energy += weight * this->get_min_max_dih_energy_ana_grad(phi_atoms.get<0>(), phi_atoms.get<1>(), phi_atoms.get<2>(), phi_atoms.get<3>(),
							phi, alpha_min_phi, alpha_max_phi,
							grad, weight * pot_weight);
					total_energy += weight * this->get_min_max_dih_energy_ana_grad(psi_atoms.get<0>(), psi_atoms.get<1>(), psi_atoms.get<2>(), psi_atoms.get<3>(),
							psi, alpha_min_psi, alpha_max_psi,
							grad, weight * pot_weight);
				}
			}

		}
	}
	return energies_map.add_energy_component(name_vector[0], total_energy);
}



}
}
}


