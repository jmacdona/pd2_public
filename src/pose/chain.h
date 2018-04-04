//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * chain.h
 *
 *  Created on: Feb 1, 2010
 *      Author: jmacdon
 */

#ifndef CHAIN_H_
#define CHAIN_H_

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>
#include "residue.h"
#include "atom.h"
#include "residue_type.h"


namespace PRODART {
namespace POSE {

class chain;
class pose;

typedef boost::shared_ptr<chain> chain_shared_ptr;
typedef boost::shared_ptr<const chain> const_chain_shared_ptr;
typedef std::vector<chain_shared_ptr> chain_shared_ptr_vector;


class chain {

	friend boost::shared_ptr<chain> new_chain();
	friend class pose;


public:

	boost::shared_ptr<pose> get_pose();
	boost::shared_ptr<const pose> get_pose() const;

	//TODO add check for this - for now always returns true
	bool isPeptide() const;
	void setIsPeptide(bool val);

	chain_shared_ptr clone() const;

	char getChainID() const;
	void setChainID(const char val);

	int length() const;
	int get_first_internal_residue_index() const;
	int get_last_internal_residue_index() const;
	int get_first_bb_atom_index() const;
	int get_last_bb_atom_index() const;

	UTILS::vector3d get_ca_pos(int internalindex) const;
	void set_ca_pos(const UTILS::vector3d pos, int internalindex) const;

	atom_shared_ptr get_ca(int internalindex);
	const_atom_shared_ptr get_ca(int internalindex) const;

	atom_shared_ptr get_bb_atom(const BBAtomType bb_at_type, int internalindex);
	const_atom_shared_ptr get_bb_atom(const BBAtomType bb_at_type, int internalindex) const;

	residue_shared_ptr get_residue(int internalindex);
	const residue_shared_ptr get_residue(int internalindex) const;

	const residue_type_vector get_sequence() const;
	//! insert a blank residue where ca-ca dist is above a cut-off
	const residue_type_vector get_sequence_with_estimated_gaps(const double cutoff = 4.1) const;

protected:

	void set_pose(boost::shared_ptr<pose> ptr);
	void set_first_internal_residue_index( const int val);
	void set_last_internal_residue_index( const int val);
	void set_first_bb_atom_index( const int val);
	void set_last_bb_atom_index( const int val);


private:

	chain();
	chain(const chain&);

	void init();

	boost::weak_ptr<pose> this_pose;

	char chainID;

	int first_internal_residue_index, last_internal_residue_index;
	int first_bb_atom_index, last_bb_atom_index;

	bool b_is_peptide;


};

boost::shared_ptr<chain> new_chain();


inline boost::shared_ptr<pose> chain::get_pose(){
	return this->this_pose.lock();

}

inline boost::shared_ptr<const pose> chain::get_pose() const{
	return this->this_pose.lock();
}



}
}
#endif /* CHAIN_H_ */
