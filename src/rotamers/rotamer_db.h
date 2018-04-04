/*
 * rotamer_db_interface.h
 *
 *  Created on: 29 Aug 2011
 *      Author: jmacdona
 */

#ifndef ROTAMER_DB_H_
#define ROTAMER_DB_H_


#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/unordered_map.hpp>


#include "pose/residue_type.h"
#include "rotamer_entry.h"
#include "pose/pose.h"


namespace PRODART {
namespace ROTAMERS {

class rotamer_db;

typedef boost::shared_ptr<rotamer_db> rotamer_db_shared_ptr;
typedef boost::shared_ptr<const rotamer_db> const_rotamer_db_shared_ptr;
typedef std::vector<rotamer_db_shared_ptr> rotamer_db_shared_ptr_vector;
typedef std::vector<const_rotamer_db_shared_ptr> const_rotamer_db_shared_ptr_vector;

typedef boost::unordered_map<int, const_rotamer_entry_shared_ptr> int_const_rotamer_entry_shared_ptr_unord_map;
typedef boost::unordered_map<int, int_const_rotamer_entry_shared_ptr_unord_map> int_int_const_rotamer_entry_shared_ptr_unord_map;
typedef boost::unordered_map<const POSE::residue_type, int_int_const_rotamer_entry_shared_ptr_unord_map> rt_int_int_const_rotamer_entry_shared_ptr_unord_map;

typedef boost::unordered_map<int, const_rotamer_entry_shared_ptr_vector> int_const_rotamer_entry_shared_ptr_vector_unord_map;
typedef boost::unordered_map<const POSE::residue_type, int_const_rotamer_entry_shared_ptr_vector_unord_map> rt_int_const_rotamer_entry_shared_ptr_vector_unord_map;

class rotamer_db {
	friend rotamer_db_shared_ptr new_rotamer_db();

private:

	rotamer_db();
	rotamer_db(const rotamer_db&);
	void clear();

	rotamer_entry_shared_ptr_vector all_rotamers;
	rt_int_int_const_rotamer_entry_shared_ptr_unord_map rt_pp_bin_r_bin_rotomers;
	rt_int_const_rotamer_entry_shared_ptr_vector_unord_map rt_pp_bin_rotomers_vec;
	bool is_bb_dep;
	double pp_bin_size_deg;
	double pp_bin_size_rad;
	double pp_ignore_val;
	bool use_ignore_val;

	void index_rotamers();

public:
	virtual bool load_data( std::istream& input ) ;
	virtual bool load_gz_data( std::istream& input );

	virtual const_rotamer_entry_shared_ptr get_rotamer_entry(const POSE::residue_type& type,
			const int combined_pp_bin,
			const int combined_r_bin) const;

	virtual const_rotamer_entry_shared_ptr get_rotamer_entry(const POSE::residue_type& type,
			const double phi,
			const double psi,
			const double chi1,
			const double chi2,
			const double chi3,
			const double chi4) const;

	virtual const const_rotamer_entry_shared_ptr_vector get_rotamer_entry_vector(const POSE::residue_type& type,
			const int combined_pp_bin) const;

	virtual const const_rotamer_entry_shared_ptr_vector get_rotamer_entry_vector(const POSE::residue_type& type,
			const double phi,
			const double psi) const;






};

rotamer_db_shared_ptr new_rotamer_db();

}
}

#endif /* ROTAMER_DB_H_ */
