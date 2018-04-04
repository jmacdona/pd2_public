//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * atom.h
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#ifndef ATOM_H_
#define ATOM_H_
#include "atom_type.h"
#include "utils/vector3d.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

namespace PRODART {
namespace POSE {

class atom;
class chain;
class pose;
class residue;
class sidechain;

typedef boost::shared_ptr<atom> atom_shared_ptr;
typedef boost::shared_ptr<const atom> const_atom_shared_ptr;
typedef std::vector<atom_shared_ptr> atom_shared_ptr_vector;
typedef std::vector<const_atom_shared_ptr> const_atom_shared_ptr_vector;

class atom{

	friend boost::shared_ptr<atom> new_atom();
	friend boost::shared_ptr<atom> new_atom(const atom_type aty, const PRODART::UTILS::vector3d& vec);
	friend class pose;
	friend class residue;
	friend class sidechain;

private:

	PRODART::UTILS::vector3d coords;
	PRODART::UTILS::vector3d backed_up_coords;

	atom_type type;

	boost::weak_ptr<pose> this_pose;
	boost::weak_ptr<chain> this_chain;
	boost::weak_ptr<residue> this_residue;
	boost::weak_ptr<sidechain> this_sidechain;

	double occupancy, b_factor;
	double charge;
	double mass;

	bool active, set, has_moved_flag;

	int seq_num;

	atom(){
		occupancy = 1.0;
		b_factor = 0;
		charge = 0;
		mass = 0;
		active = false;
		set = false;
		seq_num = 0;
		has_moved_flag = true;
	}

	atom(const atom_type aty, const PRODART::UTILS::vector3d& vec){
		occupancy = 1.0;
		b_factor = 0;
		charge = 0;
		mass = 0;
		type = aty;
		coords = vec;
		active = true;
		set = true;
		seq_num = 0;
		has_moved_flag = true;
	}

	atom(const atom&);

public:




	/*
	atom& operator=(const atom& old_atom){
		this->coords = old_atom.coords;
		this->type = old_atom.type;
		return *this;
	}
	*/


	~atom();

	atom_shared_ptr clone() const;

	boost::shared_ptr<pose> get_pose(){
		return this->this_pose.lock();
	}
	boost::shared_ptr<const pose> get_pose() const{
		return this->this_pose.lock();
	}

	boost::shared_ptr<chain> get_chain(){
		return this->this_chain.lock();
	}
	boost::shared_ptr<const chain> get_chain() const{
		return this->this_chain.lock();
	}

	boost::shared_ptr<residue> get_residue(){
		return this->this_residue.lock();
	}
	boost::shared_ptr<const residue> get_residue() const{
		return this->this_residue.lock();
	}

	boost::shared_ptr<sidechain> get_sidechain(){
		return this->this_sidechain.lock();
	}
	boost::shared_ptr<const sidechain> get_sidechain() const{
		return this->this_sidechain.lock();
	}


	PRODART::UTILS::vector3d get_coords() const{
		return coords;
	}
	double get_r(const int i ) const{
		return coords[i];
	}
	double get_x() const{
		return coords.x;
	}
	double get_y() const{
		return coords.y;
	}
	double get_z() const{
		return coords.z;
	}
	void set_coords(const PRODART::UTILS::vector3d& new_coords ){
		coords = new_coords;
	}
	void set_r(const int i, double val ) {
		coords[i] = val;
	}
	void set_x( double val) {
		coords.x = val;
	}
	void set_y( double val) {
		coords.y = val;
	}
	void set_z( double val) {
		coords.z = val;
	}

	/* UNSAFE
	PRODART::UTILS::vector3d get_backup_coords() const{
		return this->backed_up_coords;
	}
	void set_backup_coords(const PRODART::UTILS::vector3d vec) {
		backed_up_coords = vec;
	}
	*/

	bool isActive() const{
		return active;
	}

	void setActive(const bool state){
		active = state;
	}

	bool isSet() const{
		return set;
	}

	void setSet(const bool state){
		set = state;
	}

	bool isActiveAndSet() const{
		return (set && active);
	}

	atom_type get_type() const {
		return type;
	}

	void set_type(atom_type new_type){
		type = new_type;
	}

	double get_occupancy() const{
		return occupancy;
	}
	void set_occupancy(const double val){
		this->occupancy = val;
	}

	double get_b_factor() const{
		return this->b_factor;
	}
	void set_b_factor(const double val){
		this->b_factor = val;
	}

	double get_charge() const{
		return this->charge;
	}
	void set_charge(const double val){
		this->charge = val;
	}

	double get_mass() const{
		return this->mass;
	}
	void set_mass(const double val){
		this->mass = val;
	}

	/* UNSAFE
	void backup_coords(){
		this->backed_up_coords = this->coords;
	}
	void restore_backed_up_coords(){
		this->coords = this->backed_up_coords;
	}
	*/

	//! NOTE: this number is likely to be highly volatile - particularly if sidechains are being redesigned etc
	int get_seq_num() const{
		return this->seq_num;
	}

	bool has_moved() const {
		return has_moved_flag;
	}

	void set_has_moved_flag(bool flag_){
		has_moved_flag = flag_;
	}

protected:

	void set_pose(boost::shared_ptr<pose> ptr){
		this->this_pose = ptr;
	}

	void set_chain(boost::shared_ptr<chain> ptr){
		this->this_chain = ptr;
	}

	void set_residue(boost::shared_ptr<residue> ptr){
		this->this_residue = ptr;
	}

	void set_sidechain(boost::shared_ptr<sidechain> ptr){
		this->this_sidechain = ptr;
	}

	void set_seq_num(const int val){
		this->seq_num = val;
	}



};

boost::shared_ptr<atom> new_atom();
boost::shared_ptr<atom> new_atom(const atom_type aty, const PRODART::UTILS::vector3d& vec);

atom_shared_ptr find( const atom_shared_ptr_vector vec, const atom_type);

}
}




#endif /* ATOM_H_ */
