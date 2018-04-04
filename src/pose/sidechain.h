//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * sidechain.h
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#ifndef SIDECHAIN_H_
#define SIDECHAIN_H_
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <vector>
#include "atom.h"
#include <iostream>


namespace PRODART {
namespace POSE {

class residue;
class sidechain;
class residue_type;

typedef boost::shared_ptr<sidechain> sidechain_shared_ptr;
typedef boost::shared_ptr<const sidechain> const_sidechain_shared_ptr;
typedef std::vector<sidechain_shared_ptr> sidechain_shared_ptr_vector;

typedef boost::weak_ptr<sidechain> sidechain_weak_ptr;

typedef std::vector<atom_shared_ptr_vector> atom_shared_ptr_vector_vector;


class sidechain{


	friend sidechain_shared_ptr new_sidechain();
	friend class residue;
	friend class pose;

protected:
	sidechain();
	sidechain(const sidechain&);



	void init(const sidechain&);


	void init();
	boost::weak_ptr<residue> res;
	sidechain_weak_ptr _this;
	atom_shared_ptr_vector sc_atom_vec;

protected:
	atom_type_vector_vector chi_defs;
	atom_type_vector_vector chi_fwd;
	atom_type_vector_vector chi_bwd;
	atom_shared_ptr_vector_vector chi_defs_atom_cache;
	atom_shared_ptr_vector_vector chi_fwd_atom_cache;
	atom_shared_ptr_vector_vector chi_bwd_atom_cache;

	const atom_shared_ptr_vector null_atom_shared_ptr_vector;

	atom_shared_ptr_vector get_atom_vec(const atom_type_vector& vec);
	void update_atom_vec_vec(const atom_type_vector_vector& vec , atom_shared_ptr_vector_vector& ret_vec);


public:

	virtual ~sidechain(){
	}

	sidechain_shared_ptr clone() const;

	void copy_sidechain(const_sidechain_shared_ptr sc_to_copy);

	//! returns first atom of type "type"
	atom_shared_ptr get_atom(const atom_type type);
	atom_shared_ptr get_atom(const int index);

	const_atom_shared_ptr get_atom(const atom_type type) const;
	const_atom_shared_ptr get_atom(const int index) const;

	atom_shared_ptr_vector get_all_sc_atoms();
	const_atom_shared_ptr_vector get_all_sc_atoms() const;

	atom_shared_ptr add_new_atom(const PRODART::UTILS::vector3d& vec, const atom_type type);

	virtual void update_chi_atom_cache();

	virtual void set_chi_atoms(const unsigned int chi_index,
			const atom_type_vector& chi_vec,
			const atom_type_vector&  chi_f_vec,
			const atom_type_vector&  chi_b_vec){


		/*
		if (chi_index > 4 || chi_index == 0){
			std::cerr << "sidechain: ERROR: chi_index out of range" << std::endl;
			return;
		}
		*/

		if (chi_index > chi_defs.size()){
			chi_defs.resize(chi_index);
			chi_fwd.resize(chi_index);
			chi_bwd.resize(chi_index);
			chi_defs_atom_cache.resize(chi_index);
			chi_fwd_atom_cache.resize(chi_index);
			chi_bwd_atom_cache.resize(chi_index);
		}

		chi_defs[chi_index-1] = chi_vec;
		chi_fwd[chi_index-1] = chi_f_vec;
		chi_bwd[chi_index-1] = chi_b_vec;

		this->update_chi_atom_cache();

	}

	virtual double get_chi(const unsigned int chi_index) const{
		if (chi_index > chi_defs.size() || chi_index > chi_defs_atom_cache.size()){
			return 0;
		}
		if (chi_defs_atom_cache[chi_index-1].size() != 4){
			return 0;
		}
		else {
			// chi1: N-CA-CB-XG
			return UTILS::dihedral(chi_defs_atom_cache[chi_index-1][0]->get_coords(),
					chi_defs_atom_cache[chi_index-1][1]->get_coords(),
					chi_defs_atom_cache[chi_index-1][2]->get_coords(),
					chi_defs_atom_cache[chi_index-1][3]->get_coords());
		}

	};

	// residue specific
	virtual int get_chi_r(const unsigned int chi_index) const;

	virtual boost::tuple<const atom_shared_ptr_vector&, const atom_shared_ptr_vector&, const atom_shared_ptr_vector&> get_chi_details(const unsigned int chi_index){
		if (chi_index <= chi_defs.size() && chi_index <= chi_defs_atom_cache.size() && chi_index != 0){
			return boost::tuple<const atom_shared_ptr_vector&, const atom_shared_ptr_vector&, const atom_shared_ptr_vector&>(chi_defs_atom_cache[chi_index-1],
					chi_fwd_atom_cache[chi_index-1],
					chi_bwd_atom_cache[chi_index-1]);
		}
		else {
			return boost::tuple<const atom_shared_ptr_vector&, const atom_shared_ptr_vector&, const atom_shared_ptr_vector&>(null_atom_shared_ptr_vector,null_atom_shared_ptr_vector,null_atom_shared_ptr_vector);
		}
	}

	virtual bool set_chi(const unsigned int chi_index, const double rad) {
		std::cerr << "sidechain: set_chi: ERROR: use derived class sidechain_rotamer to set chi angles instead" << std::endl;
		return false;
	};
	virtual bool set_chi_inverse(const unsigned int chi_index, const double rad) {
		std::cerr << "sidechain: set_chi: ERROR: use derived class sidechain_rotamer to set chi angles instead" << std::endl;
		return false;
	};

	unsigned int get_chi_count() const{
		return chi_defs.size();
	}


	int get_atom_count() const;

	boost::shared_ptr<residue> get_residue();
	boost::shared_ptr<const residue> get_residue() const;



	void clear();



protected:
	void set_residue(const boost::shared_ptr<residue> ptr);





};

const int get_chi_r(const residue_type& rt, const unsigned int chi_index, const double chi);


sidechain_shared_ptr new_sidechain();


inline atom_shared_ptr sidechain::get_atom(const atom_type type){
	return find(sc_atom_vec, type);
}

inline atom_shared_ptr sidechain::get_atom(const int index){
	return sc_atom_vec[index];
}

inline const_atom_shared_ptr sidechain::get_atom(const atom_type type) const{
	return find(sc_atom_vec, type);
}

inline const_atom_shared_ptr sidechain::get_atom(const int index) const{
	return sc_atom_vec[index];
}

inline int sidechain::get_atom_count() const{
	return static_cast<int>(sc_atom_vec.size());
}

}
}


#endif /* SIDECHAIN_H_ */
