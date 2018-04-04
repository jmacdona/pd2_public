//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * sidechain.cpp
 *
 *  Created on: Feb 1, 2010
 *      Author: jmacdon
 */

#include "sidechain.h"
#include "residue.h"
#include "pose.h"

using PRODART::UTILS::vector3d;

namespace PRODART {
namespace POSE {



sidechain::sidechain() : null_atom_shared_ptr_vector(0){
	this->init();
}

sidechain::sidechain(const sidechain& old_sc) : null_atom_shared_ptr_vector(0){
	this->init();
}

void sidechain::init(const sidechain& old_sc){
	atom_shared_ptr_vector::const_iterator iter;
	this->init();
	for (iter = old_sc.sc_atom_vec.begin(); iter != old_sc.sc_atom_vec.end(); iter++){
		atom_shared_ptr new_sc_atom = (*iter)->clone();
		sc_atom_vec.push_back(new_sc_atom);
		new_sc_atom->set_sidechain(_this.lock());
		residue_shared_ptr currRes = this->res.lock();
		if (currRes){
			new_sc_atom->set_pose(currRes->get_pose());
			new_sc_atom->set_chain(currRes->get_chain());
			new_sc_atom->set_residue(currRes);
			new_sc_atom->set_sidechain(_this.lock());
		}

	}
	// sort out chi angle copying - done
	chi_defs = old_sc.chi_defs;
	chi_fwd = old_sc.chi_fwd;
	chi_bwd = old_sc.chi_bwd;
	this->update_chi_atom_cache();
}

void sidechain::set_residue(const residue_shared_ptr ptr){
	this->res = ptr;
	atom_shared_ptr_vector::iterator iter;
	for (iter = sc_atom_vec.begin(); iter != sc_atom_vec.end(); iter++){
		residue_shared_ptr currRes = this->res.lock();
		if (currRes){
			(*iter)->set_pose(currRes->get_pose());
			(*iter)->set_chain(currRes->get_chain());
			(*iter)->set_residue(currRes);
			(*iter)->set_sidechain(_this.lock());
		}

	}

}

sidechain_shared_ptr sidechain::clone() const{
	sidechain_shared_ptr ptr(new sidechain(*this));
	ptr->_this = ptr;
	ptr->init(*this);
	return ptr;
}

void sidechain::clear(){
	sc_atom_vec.clear();
	residue_shared_ptr rptr = res.lock();
	if (rptr){
		if (res.lock()->get_pose()){
			res.lock()->get_pose()->index();
		}
	}
}

void sidechain::copy_sidechain(const_sidechain_shared_ptr sc_to_copy){
	this->clear();
	const int num_at = sc_to_copy->get_atom_count();
	for (int i = 0; i < num_at; i++){
		const_atom_shared_ptr at_to_cp = sc_to_copy->get_atom(i);
		atom_shared_ptr new_sc_atom = at_to_cp->clone();//this->add_new_atom(at_to_cp->get_coords(),at_to_cp->get_type());
		residue_shared_ptr currRes = this->res.lock();
		if (currRes){
			new_sc_atom->set_pose(currRes->get_pose());
			new_sc_atom->set_chain(currRes->get_chain());
			new_sc_atom->set_residue(currRes);
			new_sc_atom->set_sidechain(_this.lock());
		}
		sc_atom_vec.push_back(new_sc_atom);
	}

	residue_shared_ptr rptr = res.lock();
	if (rptr){
		if (res.lock()->get_pose()){
			res.lock()->get_pose()->index();
		}
	}
}

void sidechain::init(){
	sc_atom_vec.reserve(50);
	chi_defs.clear();
	chi_fwd.clear();
	chi_bwd.clear();
	chi_defs_atom_cache.clear();
	chi_fwd_atom_cache.clear();
	chi_bwd_atom_cache.clear();
	chi_defs_atom_cache.resize(4);
	chi_fwd_atom_cache.resize(4);
	chi_bwd_atom_cache.resize(4);
	chi_defs.resize(4);
	chi_fwd.resize(4);
	chi_bwd.resize(4);
}


atom_shared_ptr sidechain::add_new_atom(const PRODART::UTILS::vector3d& vec, const atom_type type){
	atom_shared_ptr new_sc_atom = this->get_atom(type);
	if (!this->get_atom(type)){
		new_sc_atom = new_atom(type, vec);
		sc_atom_vec.push_back(new_sc_atom);
		new_sc_atom->set_sidechain(_this.lock());
		residue_shared_ptr currRes = this->res.lock();
		if (currRes){
			new_sc_atom->set_pose(currRes->get_pose());
			new_sc_atom->set_chain(currRes->get_chain());
			new_sc_atom->set_residue(currRes);
			new_sc_atom->set_sidechain(_this.lock());
		}
	}
	else {
		new_sc_atom->setActive(true);
		new_sc_atom->setSet(true);
		new_sc_atom->set_coords(vec);
	}
	return new_sc_atom;
}


boost::shared_ptr<residue> sidechain::get_residue(){
	return res.lock();
}

boost::shared_ptr<const residue> sidechain::get_residue() const{
	return res.lock();
}


atom_shared_ptr_vector sidechain::get_atom_vec(const atom_type_vector& vec){
	atom_shared_ptr_vector ret_vec;
	ret_vec.reserve(vec.size());
	atom_type_vector::const_iterator it;

	residue_shared_ptr resid = this->get_residue();

	for (it = vec.begin(); it != vec.end(); it++){
		atom_shared_ptr at = this->get_atom(*it);//resid->get_atom(*it);
		if (!at){
			at = resid->get_atom(*it);
		}
		ret_vec.push_back(at);
	}
	return ret_vec;
}

void sidechain::update_atom_vec_vec(const atom_type_vector_vector& vec , atom_shared_ptr_vector_vector& ret_vec){

	atom_type_vector_vector::const_iterator it;

	ret_vec.clear();

	for (it = vec.begin(); it != vec.end(); it++){
		ret_vec.push_back(get_atom_vec(*it));
	}

}

void sidechain::update_chi_atom_cache(){

	update_atom_vec_vec( chi_defs, chi_defs_atom_cache);
	update_atom_vec_vec( chi_fwd, chi_fwd_atom_cache);
	update_atom_vec_vec( chi_bwd, chi_bwd_atom_cache);

}

const int get_chi_r(const residue_type& rt, const unsigned int chi_index, const double chi){
	if (rt.is_equal3(residue_type("PRO"))
			&& chi_index == 1){
		if (chi >= UTILS::degrees_to_radians(0) && chi <= UTILS::degrees_to_radians(90)){
			return 1;
		}
		else {
			return 2;
		}
	}
	else if (rt.is_equal3(residue_type("TRP"))
			&& chi_index == 2){
		if (chi >= UTILS::degrees_to_radians(-180) && chi <= UTILS::degrees_to_radians(-60)){
			return 1;
		}
		else if (chi >= UTILS::degrees_to_radians(-60) && chi <= UTILS::degrees_to_radians(60)){
			return 2;
		}
		else {
			return 3;
		}
	}
	else if ((rt.is_equal3(residue_type("PHE")) || rt.is_equal3(residue_type("TYR")) || rt.is_equal3(residue_type("HIS")))
			&& chi_index == 2){
		if (chi >= UTILS::degrees_to_radians(30) && chi <= UTILS::degrees_to_radians(150)){
			return 1;
		}
		else {
			return 2;
		}
	}
	else if ((rt.is_equal3(residue_type("ASN")) || rt.is_equal3(residue_type("ASP")) || rt.is_equal3(residue_type("GLN")) || rt.is_equal3(residue_type("GLU")))
			&& (chi_index == 2 || chi_index == 3) ){
		if (chi >= UTILS::degrees_to_radians(30) && chi <= UTILS::degrees_to_radians(90)){
			return 1;
		}
		else if (chi >= UTILS::degrees_to_radians(-30) && chi <= UTILS::degrees_to_radians(30)){
			return 2;
		}
		else {
			return 3;
		}
	}
	else {
		if (chi >= UTILS::degrees_to_radians(0) && chi <= UTILS::degrees_to_radians(120)){
			return 1;
		}
		else if (chi >= UTILS::degrees_to_radians(-120) && chi <= UTILS::degrees_to_radians(0)){
			//swap of 2 and 3 is deliberate
			return 3;
		}
		else {
			return 2;
		}
	}
}

int sidechain::get_chi_r(const unsigned int chi_index) const{
	if (chi_index > chi_defs.size() || chi_index > chi_defs_atom_cache.size()){
		return 0;
	}
	if (chi_defs_atom_cache[chi_index-1].size() != 4){
		return 0;
	}
	else {
		const double chi = this->get_chi(chi_index);
		const residue_type rt = this->get_residue()->get_type();
		const int r = PRODART::POSE::get_chi_r(rt, chi_index, chi);
		return r;//get_chi_r(rt, chi_index, chi);
		/*
		if (rt.is_equal3(residue_type("PRO"))
				&& chi_index == 1){
			if (chi >= UTILS::degrees_to_radians(0) && chi <= UTILS::degrees_to_radians(90)){
				return 1;
			}
			else {
				return 2;
			}
		}
		else if (rt.is_equal3(residue_type("TRP"))
				&& chi_index == 2){
			if (chi >= UTILS::degrees_to_radians(-180) && chi <= UTILS::degrees_to_radians(-60)){
				return 1;
			}
			else if (chi >= UTILS::degrees_to_radians(-60) && chi <= UTILS::degrees_to_radians(60)){
				return 2;
			}
			else {
				return 3;
			}
		}
		else if ((rt.is_equal3(residue_type("PHE")) || rt.is_equal3(residue_type("TYR")) || rt.is_equal3(residue_type("HIS")))
				&& chi_index == 2){
			if (chi >= UTILS::degrees_to_radians(30) && chi <= UTILS::degrees_to_radians(150)){
				return 1;
			}
			else {
				return 2;
			}
		}
		else if ((rt.is_equal3(residue_type("ASN")) || rt.is_equal3(residue_type("ASP")) || rt.is_equal3(residue_type("GLN")) || rt.is_equal3(residue_type("GLU")))
				&& (chi_index == 2 || chi_index == 3) ){
			if (chi >= UTILS::degrees_to_radians(30) && chi <= UTILS::degrees_to_radians(90)){
				return 1;
			}
			else if (chi >= UTILS::degrees_to_radians(-30) && chi <= UTILS::degrees_to_radians(30)){
				return 2;
			}
			else {
				return 3;
			}
		}
		else {
			if (chi >= UTILS::degrees_to_radians(0) && chi <= UTILS::degrees_to_radians(120)){
				return 1;
			}
			else if (chi >= UTILS::degrees_to_radians(-120) && chi <= UTILS::degrees_to_radians(0)){
				//swap of 2 and 3 is deliberate
				return 3;
			}
			else {
				return 2;
			}
		}
		*/
	}
}


sidechain_shared_ptr new_sidechain(){
	sidechain_shared_ptr ptr(new sidechain);
	ptr->_this = ptr;
	return ptr;
}

atom_shared_ptr_vector sidechain::get_all_sc_atoms(){
	atom_shared_ptr_vector rtn_vec(sc_atom_vec);
	return rtn_vec;
}

const_atom_shared_ptr_vector sidechain::get_all_sc_atoms() const{
	const_atom_shared_ptr_vector rtn_vec;
	rtn_vec.reserve(sc_atom_vec.size());
	for (unsigned int i = 0; i < sc_atom_vec.size(); i++){
		rtn_vec.push_back(sc_atom_vec[i]);
	}
	return rtn_vec;
}



}
}


