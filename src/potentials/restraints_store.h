/*
 * restraints_store.h
 *
 *  Created on: 25 Jan 2011
 *      Author: jmacdona
 */

#ifndef RESTRAINTS_STORE_H_
#define RESTRAINTS_STORE_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "prodart_env/prodart_env.h"
#include <boost/tuple/tuple.hpp>

#include "pose/pose.h"



#include "pose_meta/pose_meta_interface.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {




//! singleton class
class restraints_store{

	/*
		typedef boost::tuple<std::string, char, POSE::atom_type,
							std::string, char, POSE::atom_type,
							double, double >
	 */
	template<int N>
	class rst_element{

	public:

		static const int n = N;

		std::string resids[N];
		char chains[N];
		atom_type at_types[N];

		double equil_val;
		double weight;

		double equil_val_lower, equil_val_upper;

		POSE::four_state_sec_struct secs;

		bool is_parallel;

		UTILS::vector3d coord; // for coord restraint, also for SSE axis restraints (start)
		UTILS::vector3d coord2; // for SSE axis restraints (end)

		int get_n() const{
			return n;
		}

	};

	typedef rst_element<1>  coord_rst_element;
	typedef rst_element<1>  sec_struct_rst_element;
	typedef rst_element<2>  bond_harmonic_rst_element;
	typedef rst_element<2>  sse_axis_rst_element; // NOTE: same as bond_harmonic_rst_element
	typedef rst_element<2>  contact_rst_element; // NOTE: same as bond_harmonic_rst_element
	typedef rst_element<3>  ang_harmonic_rst_element;
	typedef rst_element<4>  dih_harmonic_rst_element;

	typedef std::vector<coord_rst_element> coord_rst_vector;
	typedef std::vector<sec_struct_rst_element> sec_struct_rst_vector;
	typedef std::vector<bond_harmonic_rst_element> bond_harm_rst_vector;
	typedef std::vector<sse_axis_rst_element> sse_axis_rst_vector; // NOTE: same as bond_harm_rst_vector
	typedef std::vector<ang_harmonic_rst_element> ang_harm_rst_vector;
	typedef std::vector<dih_harmonic_rst_element> dih_harm_rst_vector;
	typedef std::vector<contact_rst_element> contact_rst_vector;


private:
	restraints_store(const restraints_store&);

protected:
	restraints_store();
	virtual ~restraints_store();

	coord_rst_vector coord_rst_vec;
	sec_struct_rst_vector secs_rst_vec;
	bond_harm_rst_vector bond_rst_vec;
	ang_harm_rst_vector ang_rst_vec;
	dih_harm_rst_vector dih_rst_vec;
	bond_harm_rst_vector strand_pair_vec;
	sse_axis_rst_vector sse_axis_rst_vec;
	contact_rst_vector contact_rst_vec;
	contact_rst_vector ca_GO_rst_vec;


public:

	static restraints_store* Instance();

	std::istream& load_restraints(std::istream& input);

	void add_coord_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
			const UTILS::vector3d pos, const double weight, const double equil_dist = 0);
	void add_sec_struct_rst(const std::string res1, const char ch1,
			const POSE::four_state_sec_struct secs, const double weight);
	void add_bond_harm_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
			const std::string res2, const char ch2, const POSE::atom_type at2,
			const double equil_val, const double weight);
	void add_sse_axis_rst(const std::string res1, const char ch1,
			const std::string res2, const char ch2,
			POSE::four_state_sec_struct sse_type, const double weight,
			const UTILS::vector3d start, const UTILS::vector3d end);
	void add_ang_harm_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
			const std::string res2, const char ch2, const POSE::atom_type at2,
			const std::string res3, const char ch3, const POSE::atom_type at3,
			const double equil_val, const double weight);
	void add_dih_harm_rst(const std::string res1, const char ch1, const POSE::atom_type at1,
			const std::string res2, const char ch2, const POSE::atom_type at2,
			const std::string res3, const char ch3, const POSE::atom_type at3,
			const std::string res4, const char ch4, const POSE::atom_type at4,
			const double equil_val, const double weight);
	void add_strand_pair_rst(const std::string res1, const char ch1,
			const std::string res2, const char ch2,
			const bool is_parallel, const double weight);
	void add_contact_rst(const std::string res1, const char ch1,
			const std::string res2, const char ch2,
			const double equil_val_lower, const double equil_val_higher,
			const double weight);

	void add_ca_GO_rst(const std::string res1, const char ch1,
			const std::string res2, const char ch2,
			const double equil_val_lower, const double equil_val_higher,
			const double weight);

	void add_rst(const std::string& rst_line, const long lineNum = -1);

	void apply_restraints(PRODART::POSE::META::pose_meta_shared_ptr pose_meta_);


};



class restraints_init_exception: public std::exception{
  virtual const char* what() const throw()
  {
    return "restraints_init_exception";
  }
};


}
}
}





#endif /* RESTRAINTS_STORE_H_ */
