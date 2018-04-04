/*
 * ca_fragment_classifier.h
 *
 *  Created on: 25 Aug 2010
 *      Author: jmacdona
 */

#ifndef CA_FRAGMENT_CLASSIFIER_H_
#define CA_FRAGMENT_CLASSIFIER_H_

#include "pose_meta/pose_meta_interface.h"
#include "pose_meta/pose_meta_defs.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include "prodart_env/prodart_env.h"
#include <exception>
#include "fragment_classifier_interface.h"
#include "pose_meta/ca_pose_meta.h"
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>

namespace PRODART {
namespace POSE {
namespace META {
namespace FRAG_CLASS {

typedef std::map<int,int, std::less< int > > int_int_map;
typedef std::map< int, double, std::less< int > > int_double_map;

typedef boost::unordered_map<int,int > int_int_umap;
typedef boost::unordered_map< int, double> int_double_umap;

class t_frag_count_element{

public:
	t_frag_count_element(void);

	~t_frag_count_element(void);


	int_double_map ss_bin_count;
	int_double_map frag_number_count;
	double total_count;


};

//typedef std::map< int, t_frag_count_element, std::less< int > > int_frag_count_ele_map;
typedef boost::unordered_map< int, t_frag_count_element> int_frag_count_ele_umap;



fragment_classifier_shared_ptr new_ca_fragment_classifier();

//! not really a singleton class but has one main instance
class ca_fragment_classifier : public fragment_classifier_interface {

	friend fragment_classifier_shared_ptr new_ca_fragment_classifier();


private:

	ca_fragment_classifier(const ca_fragment_classifier&);


protected:

	void fillGapsInIBBB_map();
	double IBBB_key_distance(const int key1, const int key2) const;
	int findClosestKey(const int search_key) const;

	static const double four_CA_bin_size;


	int_frag_count_ele_umap IBBB_SS_count;
	int_int_umap frag_number_to_SStype_map;
	int_int_map IBBB_to_frag_number_map;



public:

	ca_fragment_classifier();
	virtual ~ca_fragment_classifier(){
	}

	bool classify_fragments (const pose_shared_ptr _pose,
			const PRODART::POSE::META::ca_pose_meta_shared_ptr _pose_meta) const;

	static fragment_classifier_shared_ptr MainInstance();

	std::istream& loadData( std::istream& input);


};










}
}
}
}
#endif /* CA_FRAGMENT_CLASSIFIER_H_ */

