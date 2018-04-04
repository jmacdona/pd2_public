/*
 * backbone_builder.h
 *
 *  Created on: 15 Sep 2010
 *      Author: jmacdona
 */

#ifndef BACKBONE_BUILDER_H_
#define BACKBONE_BUILDER_H_
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <iostream>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include "prodart_env/prodart_env.h"
#include <exception>
#include "pose_meta/ca_pose_meta.h"
#include "pose_meta/frag_classify/fragment_classifier_interface.h"
#include "pose_meta/frag_classify/ca_fragment_classifier.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include "pose_utils/pose_utils.h"
#include <boost/unordered_map.hpp>
#include <utility>


namespace PRODART {
namespace POSE {
namespace BB_BUILDER {

typedef std::vector< double > double_vector;
typedef std::vector< int > int_vector;
typedef std::map< unsigned long, unsigned long > ulong_ulong_map;
typedef std::map< int, unsigned long > int_ulong_map;
typedef std::map< int, PRODART::UTILS::vector3d > int_vector3d_map;
typedef std::vector< std::string >  split_vector_type;


class backbone_builder_dih_element {

public:

	backbone_builder_dih_element(){
		count = 0;
	}

	void calcAverages();

	PRODART::UTILS::vector3d C_local, O_local, N_local, CB_local;

	double count;

private:

};

typedef std::map< int, backbone_builder_dih_element > int_bb_dih_ele_map;

class backbone_builder_key_element {

public:

	backbone_builder_key_element(){
		dih_ele_map.clear();
		total_count = 0;

	}

	void add_count_to_db(
			const PRODART::UTILS::vector3d C_local,
			const PRODART::UTILS::vector3d O_local,
			const PRODART::UTILS::vector3d N_local,
			const PRODART::UTILS::vector3d CB_local,
			const int first_dih_bin,
			const int last_dih_bin);

	void get_local_vecs(
			PRODART::UTILS::vector3d &C_local,
			PRODART::UTILS::vector3d &O_local,
			PRODART::UTILS::vector3d &N_local,
			PRODART::UTILS::vector3d &CB_local,
			const int first_dih_bin,
			const int last_dih_bin) const;

	void calcAverages();

	int discardBadGeom();

	int getClosestKey(const int) const;

	double getDihKeyDist(const int key1, const int key2) const;

	unsigned long total_count;

	int_bb_dih_ele_map dih_ele_map;

private:


};


//class backbone_builder_dih_element;

//class backbone_builder_key_element;
typedef std::map< int, backbone_builder_key_element > int_bb_ele_map;

class backbone_builder;

typedef boost::shared_ptr<const backbone_builder> const_backbone_builder_shared_ptr;


typedef boost::unordered_map< unsigned long, unsigned long > ulong_ulong_umap;
typedef boost::unordered_map< int, unsigned long > int_ulong_umap;
typedef boost::unordered_map< int, PRODART::UTILS::vector3d > int_vector3d_umap;
typedef boost::unordered_map< int, backbone_builder_dih_element > int_bb_dih_ele_umap;
typedef boost::unordered_map< int, backbone_builder_key_element > int_bb_ele_umap;


class backbone_builder {




	friend std::ostream &operator<<( std::ostream&,  backbone_builder & );
	//friend std::istream &operator>>( std::istream&, backbone_builder & );




private:

	backbone_builder(const backbone_builder&);

	static void Init();

	//MTRand *randomNumGen;
	mutable MTRand::MTRand_shared_ptr rand_gen;


	// ************
	// IBBB data store:
	//
	int_ulong_umap d_m1_p1_count_map, d_0_p2_count_map,
		d_m1_p2_count_map;

	int_ulong_umap int_overall_bin_count, CB_int_overall_bin_count;

	int_vector3d_umap C_vec_map, O_vec_map, N_vec_map, CB_vec_map;


	int_ulong_umap PRO_d_m1_p1_count_map, PRO_d_0_p2_count_map,
		PRO_d_m1_p2_count_map;

	int_ulong_umap PRO_int_overall_bin_count;

	int_vector3d_umap PRO_C_vec_map, PRO_O_vec_map, PRO_N_vec_map, PRO_CB_vec_map;
	//******************

	int_bb_ele_umap non_PRO_data_map;
	int_bb_ele_umap PRO_data_map;

	inline int getKey(double sign,
			double d_m1_p2,
			double d_0_p2,
			double d_m1_p1,
			double first_dih,
			double last_dih) const;
	inline int getKey(int IBBB_key,
			double first_dih,
			double last_dih) const;
	inline void decomposeKey(const int key,
			int &sign,
			int &d_m1_p2_bin,
			int &d_0_p2_bin,
			int &d_m1_p1_bin,
			int &first_dih_bin,
			int &last_dih_bin) const;

	double key_distance(const int key1, const int key2) const;
	double IBBB_key_distance(const int key1, const int key2) const;

	bool validate_residues(
			PRODART::POSE::const_residue_shared_ptr, PRODART::POSE::const_residue_shared_ptr,
			PRODART::POSE::const_residue_shared_ptr, PRODART::POSE::const_residue_shared_ptr,
			PRODART::POSE::const_residue_shared_ptr, PRODART::POSE::const_residue_shared_ptr);

	double bin_size;
	const double dih_bin_size;



	const int IBBB_key_multiplier, d_m1_p2_bin_multiplier, d_0_p2_bin_multiplier, d_m1_p1_bin_multiplier, first_dih_bin_multiplier, last_dih_bin_multiplier;

	const double dih_dist_weight;

	const PRODART::POSE::atom_type ref_CA, ref_CB, ref_C, ref_N, ref_O, ref_H;


	//const int num_dih_bins;
protected:
	backbone_builder();
	//virtual ~backbone_builder();


public:

	static const_backbone_builder_shared_ptr Instance();

	virtual ~backbone_builder();


	std::istream& load_bb_db( std::istream& input);
	std::istream& load_bb_IBBB_db( std::istream& input);


	//void addToDB(PRODART::POSE::META::ca_pose_meta_shared_ptr pose_meta_);

	bool buildBackbone(PRODART::POSE::META::ca_pose_meta_shared_ptr pose_meta_) const;

	bool buildBackbone_single(PRODART::POSE::META::ca_pose_meta_shared_ptr pose_meta_,
			PRODART::POSE::META::frag4_element& ele,
			const bool change_res1_CB = true,
			const bool withNoise = false,
			const double noise_level = 20) const;

	//bool addPseudoAtoms2Ends(PRODART::PROT::Pdb&);
	//bool deleteEndAtoms(PRODART::PROT::Pdb&);

	//void test(Pdb&);

	//void outputBins();

	void calcAverages();

	void discardBadGeom();
	void discardRareGeom(const unsigned long  count_cuttoff, const unsigned long  PRO_count_cutoff);
	//void outputAverageCoors();

	int findClosestKey(const int search_key) const;
	int PRO_findClosestKey(const int search_key) const;

	int findClosestKey_with_noise(const int search_key, const double dist) const;
	int PRO_findClosestKey_with_noise(const int search_key, const double dist) const;



	int IBBBfindClosestKey(const int) const;
	int PRO_IBBBfindClosestKey(const int) const;

	int IBBBfindClosestKey_with_noise(const int search_key, const double noise_dist) const;
	int PRO_IBBBfindClosestKey_with_noise(const int search_key, const double noise_dist) const;

	//IBBB_collect_data initial_backbone_builder;


	void set_rand_num_gen(MTRand::MTRand_shared_ptr ptr) const{
		rand_gen = ptr;
	}


};






//std::istream &operator>>( std::istream&, backbone_builder & );



}
}
}



#endif /* BACKBONE_BUILDER_H_ */
