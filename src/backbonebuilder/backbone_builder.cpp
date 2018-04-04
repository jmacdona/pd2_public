/*
 * backbone_builder.cpp
 *
 *  Created on: 15 Sep 2010
 *      Author: jmacdona
 */

#include "backbone_builder.h"


using namespace PRODART::UTILS;
using namespace PRODART;
using namespace PRODART::POSE;
using namespace PRODART::POSE_UTILS;
using namespace std;
using namespace boost;
using namespace PRODART::POSE::META;

namespace PRODART {
namespace POSE {
namespace BB_BUILDER {

typedef boost::shared_ptr<backbone_builder> backbone_builder_shared_ptr;
//const int NO_FRAG = 99999;



namespace {
backbone_builder_shared_ptr instance;
boost::once_flag once_flag = BOOST_ONCE_INIT;
}

void backbone_builder::Init(){
	if (!instance){
		instance = backbone_builder_shared_ptr(new backbone_builder());
	}
}

const_backbone_builder_shared_ptr backbone_builder::Instance(){
	/*
	if (!instance){
		instance = backbone_builder_shared_ptr(new backbone_builder());
	}
	*/
	boost::call_once(&backbone_builder::Init, once_flag);
	return instance;
}




void backbone_builder_dih_element::calcAverages(){

	this->C_local = this->C_local / this->count;
	this->O_local = this->O_local / this->count;
	this->N_local = this->N_local / this->count;
	this->CB_local = this->CB_local / this->count;

}




void backbone_builder_key_element::add_count_to_db(
		const PRODART::UTILS::vector3d C_local,
		const PRODART::UTILS::vector3d O_local,
		const PRODART::UTILS::vector3d N_local,
		const PRODART::UTILS::vector3d CB_local,
		const int first_dih_bin,
		const int last_dih_bin){

	const int overall_dih_key = (first_dih_bin * 100) + last_dih_bin;

	this->dih_ele_map[overall_dih_key].count++;
	this->dih_ele_map[overall_dih_key].C_local += C_local;
	this->dih_ele_map[overall_dih_key].O_local += O_local;
	this->dih_ele_map[overall_dih_key].N_local += N_local;
	this->dih_ele_map[overall_dih_key].CB_local += CB_local;

	this->total_count++;

}

int backbone_builder_key_element::discardBadGeom(){

	int_bb_dih_ele_map::iterator iter;
	int_ulong_map erase_list;

	for (iter = dih_ele_map.begin(); iter != dih_ele_map.end(); iter++){
		PRODART::UTILS::vector3d CA_vec(0,0,0);
		PRODART::UTILS::vector3d C_vec = iter->second.C_local;
		PRODART::UTILS::vector3d O_vec = iter->second.O_local;
		PRODART::UTILS::vector3d N_vec = iter->second.N_local;

		double CA_C_dist = (C_vec - CA_vec).mod();
		double C_O_dist = (C_vec - O_vec).mod();
		double C_N_dist = (C_vec - N_vec).mod();

		double CA_C_err = CA_C_dist - 1.51;
		double C_O_err = C_O_dist - 1.24;
		double C_N_err = C_N_dist - 1.33;

		if (!(fabs(CA_C_err) < 0.2 &&
						fabs(C_O_err) < 0.2 &&
						fabs(C_N_err) < 0.2 )){


			total_count -= static_cast<unsigned long >(iter->second.count);
			erase_list[iter->first] = 0;


		}

	}

	int_ulong_map::iterator iu_iter;

	for (iu_iter = erase_list.begin(); iu_iter != erase_list.end(); iu_iter++){
		dih_ele_map.erase(iu_iter->first);
	}


	return static_cast<int>(dih_ele_map.size());
}

void backbone_builder_key_element::calcAverages(){

	int_bb_dih_ele_map::iterator iter;

	for (iter = dih_ele_map.begin(); iter != dih_ele_map.end(); iter++){
		iter->second.calcAverages();
	}


}

int backbone_builder_key_element::getClosestKey(const int search_key) const{
	int_bb_dih_ele_map::const_iterator iter;
	double best_dist = DBL_MAX;
	int best_key = 0;
	for (iter = dih_ele_map.begin(); iter != dih_ele_map.end(); iter++){
		const double dist = getDihKeyDist(search_key, iter->first);
		if(dist == 0){
			return iter->first;
		}
		else if (dist < best_dist ){
            best_dist = dist;
            best_key = iter->first;
		}
	}
	return best_key;

}

double backbone_builder_key_element::getDihKeyDist(const int key1, const int key2) const{

	const int first_dih_bin1 = (key1 / 100);
	const int last_dih_bin1 = key1 - (100*first_dih_bin1);

	const int first_dih_bin2 = (key2 / 100);
	const int last_dih_bin2 = key2 - (100*first_dih_bin2);

	const double dist = sqrt(static_cast<double>((first_dih_bin1 - first_dih_bin2) * (first_dih_bin1 - first_dih_bin2))
						+ static_cast<double>((last_dih_bin1 - last_dih_bin2) * (last_dih_bin1 - last_dih_bin2)));

	return dist;
}

void backbone_builder_key_element::get_local_vecs(
		PRODART::UTILS::vector3d &C_local,
		PRODART::UTILS::vector3d &O_local,
		PRODART::UTILS::vector3d &N_local,
		PRODART::UTILS::vector3d &CB_local,
		const int first_dih_bin,
		const int last_dih_bin) const{

	const int overall_dih_key = (first_dih_bin * 100) + last_dih_bin;

	const int used_key = this->getClosestKey(overall_dih_key);

	int_bb_dih_ele_map::const_iterator iter = this->dih_ele_map.find(used_key) ;

	if (iter != dih_ele_map.end()){
		C_local = iter->second.C_local;
		O_local = iter->second.O_local;
		N_local = iter->second.N_local;
		CB_local = iter->second.CB_local;
	}

}






backbone_builder::backbone_builder() : bin_size(0.3),
										dih_bin_size((2.0 * PI) / static_cast<double>(5)),
										IBBB_key_multiplier(100),
										d_m1_p2_bin_multiplier(1000000),
										d_0_p2_bin_multiplier(10000),
										d_m1_p1_bin_multiplier(100),
										first_dih_bin_multiplier(10),
										last_dih_bin_multiplier(1),
										dih_dist_weight(1.5),
										ref_CA("CA"), ref_CB("CB"), ref_C("C"), ref_N("N"), ref_O("O"), ref_H("H"){

	const string ibbb_db_path = PRODART::ENV::get_option_value<string>("database:path:bb_builder_IBBB");
	std::ifstream input_ibbb_db(ibbb_db_path.c_str(), ios::in);
	cout << "backbone_builder: loading IBBB db" << endl;
	this->load_bb_IBBB_db(input_ibbb_db);

	const string db_path = PRODART::ENV::get_option_value<string>("database:path:backbone_builder");
	std::ifstream inputdb(db_path.c_str(), ios::in);
	cout << "backbone_builder: loading db" << endl;
	this->load_bb_db(inputdb);

	this->rand_gen = PRODART::ENV::prodart_env::Instance()->get_random_num_gen();

}


backbone_builder::~backbone_builder(){

	//cout << "backbone_builder: being destroyed" << endl;


}

inline int backbone_builder::getKey(const double sign,
		const double d_m1_p2,
		const double d_0_p2,
		const double d_m1_p1,
		const double first_dih,
		const double last_dih) const{

	const int d_m1_p1_bin = static_cast<int>( d_m1_p1 / bin_size );
	const int d_0_p2_bin = static_cast<int>( d_0_p2 / bin_size );
	const int d_m1_p2_bin = static_cast<int>( d_m1_p2 / bin_size );
	const int first_dih_bin = static_cast<int>( (PI + first_dih) / dih_bin_size );
	const int last_dih_bin = static_cast<int>( (PI + last_dih) / dih_bin_size );

	const int overall_key = static_cast<int>(sign) * (d_m1_p2_bin_multiplier * abs(d_m1_p2_bin)
					+ d_0_p2_bin_multiplier * d_0_p2_bin
					+ d_m1_p1_bin_multiplier *  d_m1_p1_bin
					+ first_dih_bin_multiplier * first_dih_bin
					+ last_dih_bin_multiplier * last_dih_bin);


	/*
	cout <<  overall_key << "\t"
		 << d_m1_p2_bin << "\t"
		 <<	d_0_p2_bin << "\t"
		 <<	d_m1_p1_bin << "\t"
		 <<	first_dih_bin << "\t"
		 <<	last_dih_bin << "\t"
		 << endl;
	*/


	return overall_key;

}

inline int backbone_builder::getKey(int IBBB_key,
		double first_dih,
		double last_dih) const{

	const int sign = IBBB_key >= 0 ? 1 : -1;
	const int first_dih_bin = static_cast<int>( (PI + first_dih) / dih_bin_size );
	const int last_dih_bin = static_cast<int>( (PI + last_dih) / dih_bin_size );

	const int overall_key = static_cast<int>(sign) * (this->IBBB_key_multiplier * abs(IBBB_key)
					+ first_dih_bin_multiplier * first_dih_bin
					+ last_dih_bin_multiplier * last_dih_bin);


	/*
	cout <<  overall_key << "\t"
		 << d_m1_p2_bin << "\t"
		 <<	d_0_p2_bin << "\t"
		 <<	d_m1_p1_bin << "\t"
		 <<	first_dih_bin << "\t"
		 <<	last_dih_bin << "\t"
		 << endl;
	*/


	return overall_key;
}

inline void backbone_builder::decomposeKey(const int key,
		int &sign,
		int &d_m1_p2_bin,
		int &d_0_p2_bin,
		int &d_m1_p1_bin,
		int &first_dih_bin,
		int &last_dih_bin) const{

	const int abs_key = abs(key);

	sign = key / abs_key;


	d_m1_p2_bin = (abs_key / d_m1_p2_bin_multiplier);
	d_0_p2_bin = (abs_key - abs(d_m1_p2_bin * d_m1_p2_bin_multiplier))
								/ d_0_p2_bin_multiplier;
	d_m1_p1_bin = (abs_key - abs(d_m1_p2_bin * d_m1_p2_bin_multiplier)
								- (d_0_p2_bin * d_0_p2_bin_multiplier)) / d_m1_p1_bin_multiplier;

	first_dih_bin = (abs_key - abs(d_m1_p2_bin * d_m1_p2_bin_multiplier)
			- (d_0_p2_bin * d_0_p2_bin_multiplier)
			- (d_m1_p1_bin * d_m1_p1_bin_multiplier)) / first_dih_bin_multiplier;

	last_dih_bin = (abs_key - abs(d_m1_p2_bin * d_m1_p2_bin_multiplier)
			- (d_0_p2_bin * d_0_p2_bin_multiplier)
			- (d_m1_p1_bin * d_m1_p1_bin_multiplier)
			- (first_dih_bin * first_dih_bin_multiplier)) / last_dih_bin_multiplier;

}


/*
void backbone_builder::addToDB( PRODART::POSE::META::ca_pose_meta_shared_ptr pose_meta_ ){


	const string dash("-");
	const string space(" ");

	protein.addCBtoGLY();
	protein.index();

	string resTyStr_0;


	string extChID;
	string ss_type;

	int resCount = protein.getResidueCount();


	for (int i = 2; i < resCount - 3 ; i++){
		//Residue *res_i_minus1, *res_i, *res_iplus1, *res_iplus2;

		Residue *const res_i_minus2 = protein.getResidue(i - 2);
		Residue *const res_i_minus1 = protein.getResidue(i - 1);
		Residue *const res_i = protein.getResidue(i);
		Residue *const res_iplus1 = protein.getResidue(i+1);
		Residue *const res_iplus2 = protein.getResidue(i+2);
		Residue *const res_iplus3 = protein.getResidue(i+3);



		if (this->validate_residues(
				res_i_minus2,
				res_i_minus1,
				res_i, res_iplus1,
				res_iplus2,
				res_iplus3)){
			myVector CA_vecs[4], C_vecs[4], O_vecs[4], N_vecs[4], CB_vecs[4];

			ResidueType resTy_p1 = res_iplus1->getResidueType(0);
			//ResidueType resTy_0 = res_i->getResidueType(0);

			for ( int a = 0; a < 4; a++  ){
				CA_vecs[a] = protein.getAtom( i - 1 + a, CA)->getVec();
				C_vecs[a] = protein.getAtom( i - 1 + a, C)->getVec();
				O_vecs[a] = protein.getAtom( i - 1 + a, O)->getVec();
				N_vecs[a] = protein.getAtom( i - 1 + a, N)->getVec();
			}

			const myVector CA_i_m2 = protein.getAtom( i - 2, CA)->getVec();
			const myVector CA_i_p3 = protein.getAtom( i + 3, CA)->getVec();

			//if (resTy_0 != GLY){
			CB_vecs[1] = protein.getAtom( i, CB)->getVec();
			//}


			const myVector v_i_minus1 = CA_vecs[1] - CA_vecs[0];
			const myVector v_i = CA_vecs[2] - CA_vecs[1];
			const myVector v_i_plus1 = CA_vecs[3] - CA_vecs[2];

			const myVector v_p = CA_vecs[3] - CA_vecs[1];

			double chi = (v_i_minus1 * v_i).dot(v_i_plus1);

			chi = chi / fabs(chi);

			const double d_m1_p1 = (CA_vecs[0] - CA_vecs[2]).mod();
			const double d_0_p2 = (CA_vecs[1] - CA_vecs[3]).mod();
			const double d_m1_p2 = chi * (CA_vecs[0] - CA_vecs[3]).mod();

			const double dih_first = dihedral(CA_i_m2, CA_vecs[0],CA_vecs[1],CA_vecs[2] );
			const double dih_last = dihedral(CA_vecs[1],CA_vecs[2],CA_vecs[3],CA_i_p3 );

			const int first_dih_bin = static_cast<int>( (PI + dih_first) / dih_bin_size );
			const int last_dih_bin = static_cast<int>( (PI + dih_last) / dih_bin_size );



			const int overall_key = this->getKey(chi,
					d_m1_p2,
					d_0_p2,
					d_m1_p1,
					dih_first,
					dih_last);



			const myVector x_unit = (v_i * v_p) / (v_i * v_p).mod();
			const myVector y_unit = (v_p * x_unit) / (v_p * x_unit).mod();
			const myVector z_unit = ( x_unit * y_unit ) / ( x_unit * y_unit ).mod();

			const myVector C_i_vec = C_vecs[1] - CA_vecs[1];
			const myVector CB_i_vec = CB_vecs[1] - CA_vecs[1];
			const myVector O_i_vec = O_vecs[1] - CA_vecs[1];
			const myVector N_i_p1_vec = N_vecs[2] - CA_vecs[1];

			myVector C_local, O_local, N_local, CB_local;

			C_local.x = x_unit.dot(C_i_vec);
			C_local.y = y_unit.dot(C_i_vec);
			C_local.z = z_unit.dot(C_i_vec);

			CB_local.x = x_unit.dot(CB_i_vec);
			CB_local.y = y_unit.dot(CB_i_vec);
			CB_local.z = z_unit.dot(CB_i_vec);

			O_local.x = x_unit.dot(O_i_vec);
			O_local.y = y_unit.dot(O_i_vec);
			O_local.z = z_unit.dot(O_i_vec);

			N_local.x = x_unit.dot(N_i_p1_vec);
			N_local.y = y_unit.dot(N_i_p1_vec);
			N_local.z = z_unit.dot(N_i_p1_vec);



			if (resTy_p1 != PRO){


				non_PRO_data_map[overall_key].add_count_to_db(C_local,
						O_local,
						N_local,
						CB_local,
						first_dih_bin,
						last_dih_bin);


			}
			else {

				PRO_data_map[overall_key].add_count_to_db(C_local,
						O_local,
						N_local,
						CB_local,
						first_dih_bin,
						last_dih_bin);

			}


		}


	}

}
*/

void backbone_builder::calcAverages(){


	int_bb_ele_umap::iterator iter;

	for ( iter = non_PRO_data_map.begin(); iter != non_PRO_data_map.end(); iter++){
		iter->second.calcAverages();
	}

	for ( iter = PRO_data_map.begin(); iter != PRO_data_map.end(); iter++){
		iter->second.calcAverages();
	}


	/*
	int_ulong_map::const_iterator iter;

	for ( iter = int_overall_bin_count.begin(); iter != int_overall_bin_count.end(); iter++ ){
		C_vec_map[iter->first] = C_vec_map[iter->first] / static_cast<double>(iter->second);
		O_vec_map[iter->first] = O_vec_map[iter->first] / static_cast<double>(iter->second);
		N_vec_map[iter->first] = N_vec_map[iter->first] / static_cast<double>(iter->second);
		CB_vec_map[iter->first] = CB_vec_map[iter->first] / static_cast<double>(iter->second);
	}



	for ( iter = PRO_int_overall_bin_count.begin(); iter != PRO_int_overall_bin_count.end(); iter++ ){
		PRO_C_vec_map[iter->first] = PRO_C_vec_map[iter->first] / static_cast<double>(iter->second);
		PRO_O_vec_map[iter->first] = PRO_O_vec_map[iter->first] / static_cast<double>(iter->second);
		PRO_N_vec_map[iter->first] = PRO_N_vec_map[iter->first] / static_cast<double>(iter->second);
		PRO_CB_vec_map[iter->first] = PRO_CB_vec_map[iter->first] / static_cast<double>(iter->second);
	}
	*/

}

void backbone_builder::discardRareGeom(const unsigned long  count_cuttoff, const unsigned long  PRO_count_cutoff){

	int_bb_ele_umap::iterator iter;
	int_ulong_umap reg_erase_list;
	int_ulong_umap PRO_erase_list;

	for (iter = non_PRO_data_map.begin(); iter != non_PRO_data_map.end(); iter++){
		//const int count =  iter->second.discardBadGeom();
		if (iter->second.total_count <=  count_cuttoff){
			reg_erase_list[iter->first] = 0;
		}

	}

	for (iter = PRO_data_map.begin(); iter != PRO_data_map.end(); iter++){
		//const int count =  iter->second.discardBadGeom();
		if (iter->second.total_count <=  PRO_count_cutoff){
			PRO_erase_list[iter->first] = 0;
		}

	}


	int_ulong_umap::iterator iu_iter;

	for (iu_iter = reg_erase_list.begin(); iu_iter != reg_erase_list.end(); iu_iter++){
		non_PRO_data_map.erase(iu_iter->first);
	}

	for (iu_iter = PRO_erase_list.begin(); iu_iter != PRO_erase_list.end(); iu_iter++){
		PRO_data_map.erase(iu_iter->first);
	}


	std::cerr << "non PRO discared:\t" <<  reg_erase_list.size() << "\n";
	std::cerr << "PRO discared:\t" <<  PRO_erase_list.size() << "\n";


}

void backbone_builder::discardBadGeom(){

	int_bb_ele_umap::iterator iter;
	int_ulong_umap reg_erase_list;
	int_ulong_umap PRO_erase_list;

	for (iter = non_PRO_data_map.begin(); iter != non_PRO_data_map.end(); iter++){
		const int count =  iter->second.discardBadGeom();
		if (count == 0 ){
			reg_erase_list[iter->first] = 0;
		}

	}

	for (iter = PRO_data_map.begin(); iter != PRO_data_map.end(); iter++){
		const int count =  iter->second.discardBadGeom();
		if (count == 0 ){
			PRO_erase_list[iter->first] = 0;
		}

	}


	int_ulong_umap::iterator iu_iter;

	for (iu_iter = reg_erase_list.begin(); iu_iter != reg_erase_list.end(); iu_iter++){
		non_PRO_data_map.erase(iu_iter->first);
	}

	for (iu_iter = PRO_erase_list.begin(); iu_iter != PRO_erase_list.end(); iu_iter++){
		PRO_data_map.erase(iu_iter->first);
	}




}






bool backbone_builder::buildBackbone(PRODART::POSE::META::ca_pose_meta_shared_ptr pose_meta_) const{

	const pose_shared_ptr pose_ = pose_meta_->get_pose();


	//int prev_frag_num = NO_FRAG;
	//int frag_num = NO_FRAG;
	//int next_frag_num = NO_FRAG;


	frag4_vector& fragments =  pose_meta_->get_fragments();

	//frag4_vector::const_iterator iter;
	for (unsigned int i = 0; i < fragments.size(); i++){
		//frag_num = fragments[i].frag_type_num;
		frag4_element& ele = fragments[i];

		/*
		if (i+1 < fragments.size()){
			next_frag_num = fragments[i+1].frag_type_num;
		}
		if (fragments[i].frag_pos == fragNTERM || fragments[i].frag_pos == fragNandCTERM){
			prev_frag_num = NO_FRAG;
		}
		if (fragments[i].frag_pos == fragCTERM || fragments[i].frag_pos == fragNandCTERM){
			next_frag_num = NO_FRAG;
		}
		*/





		this->buildBackbone_single(pose_meta_, ele);



		//prev_frag_num = frag_num;



	}

	// TODO use better method
	quick_add_HN(pose_, true);

	return true;
}

bool backbone_builder::buildBackbone_single(PRODART::POSE::META::ca_pose_meta_shared_ptr pose_meta_,
		PRODART::POSE::META::frag4_element& ele,
		const bool changeCB,
		const bool withNoise,
		const double noise_level) const{

	const pose_shared_ptr pose_ = pose_meta_->get_pose();



	const int resNum_1 = ele.CA_1_residue_number;
	const int resNum_m2 = ele.CA_1_residue_number - 2;
	const int resNum_p3 = ele.CA_1_residue_number + 3;

	atom_shared_ptr ca_m2 = atom_shared_ptr(), ca_p3 = atom_shared_ptr();

	if (resNum_m2 >= 0
			&& resNum_p3 < pose_->get_residue_count()  ){

		ca_m2 = pose_->get_bb_atom(CA, resNum_m2);
		ca_p3 = pose_->get_bb_atom(CA,resNum_p3);

		if (!(ca_m2->isActive() && ca_p3->isActive())){
			ca_m2 = atom_shared_ptr();
			ca_p3 = atom_shared_ptr();
		}
		else if (!(ca_m2->get_chain() == ele.CA_atoms[1]->get_chain()
				&& ca_p3->get_chain() == ele.CA_atoms[1]->get_chain())){
			ca_m2 = atom_shared_ptr();
			ca_p3 = atom_shared_ptr();
		}

	}


	if (ca_m2 && ca_p3){
		const double dih_first = dihedral(ca_m2->get_coords(),
				ele.CA_atoms[0]->get_coords(),
				ele.CA_atoms[1]->get_coords(),
				ele.CA_atoms[2]->get_coords() );
		const double dih_last = dihedral(ele.CA_atoms[1]->get_coords(),
				ele.CA_atoms[2]->get_coords(),
				ele.CA_atoms[3]->get_coords(),
				ca_p3->get_coords() );

		const int first_dih_bin = static_cast<int>( (PI + dih_first) / dih_bin_size );
		const int last_dih_bin = static_cast<int>( (PI + dih_last) / dih_bin_size );

		const int overall_key = this->getKey(ele.IBBB_key,
				dih_first,
				dih_last);


		//const vector3d v_i_minus1 = ele.CA_atoms[1]->get_coords() - ele.CA_atoms[0]->get_coords();
		const vector3d v_i = ele.CA_atoms[2]->get_coords() - ele.CA_atoms[1]->get_coords();
		//const vector3d v_i_plus1 = ele.CA_atoms[3]->get_coords() - ele.CA_atoms[2]->get_coords();
		const vector3d v_p = ele.CA_atoms[3]->get_coords() - ele.CA_atoms[1]->get_coords();

		const vector3d x_unit = (v_i * v_p) / (v_i * v_p).mod();
		const vector3d y_unit = (v_p * x_unit) / (v_p * x_unit).mod();
		const vector3d z_unit = ( x_unit * y_unit ) / ( x_unit * y_unit ).mod();

		if (!(pose_->get_residue(resNum_1)->get_type() == residue_type("PRO"))){

			int used_key;
			if (withNoise){
				used_key = this->findClosestKey_with_noise(overall_key, noise_level);
			}
			else {
				used_key = findClosestKey(overall_key);
			}

			vector3d C_local;
			vector3d CB_local;
			vector3d O_local;
			vector3d N_local;

			non_PRO_data_map.find(used_key)->second.get_local_vecs(
					C_local,
					O_local,
					N_local,
					CB_local,
					first_dih_bin,
					last_dih_bin);

			const vector3d C_vec = ele.CA_atoms[1]->get_coords() + (C_local.x * x_unit
					+ C_local.y * y_unit
					+ C_local.z * z_unit);
			const vector3d CB_vec = ele.CA_atoms[1]->get_coords() + (CB_local.x * x_unit
					+ CB_local.y * y_unit
					+ CB_local.z * z_unit);
			const vector3d O_vec = ele.CA_atoms[1]->get_coords() + (O_local.x * x_unit
					+ O_local.y * y_unit
					+ O_local.z * z_unit);
			const vector3d N_vec = ele.CA_atoms[1]->get_coords() + (N_local.x * x_unit
					+ N_local.y * y_unit
					+ N_local.z * z_unit);

			if (!pose_->get_bb_atom(POSE::C, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "C" << endl;
			}
			pose_->add_new_atom(C_vec, ref_C, resNum_1);
			if (!(pose_->get_residue(resNum_1)->get_type() == residue_type("GLY")) && changeCB){
				if (!pose_->get_bb_atom(POSE::CB, resNum_1)->isSet()){
					cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
							<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
							<< "CB" << endl;
				}
				pose_->add_new_atom(CB_vec, ref_CB, resNum_1);
			}
			if (!pose_->get_bb_atom(POSE::O, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "O" << endl;
			}

			pose_->add_new_atom(O_vec, ref_O, resNum_1);
			if (!pose_->get_bb_atom(POSE::N, resNum_1+1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1+1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1+1)->get_chain()->getChainID() << "\t"
						<< "N" << endl;
			}
			pose_->add_new_atom(N_vec, ref_N, resNum_1+1);


			/*
			res_i->newAtom(C, C_vec.x, C_vec.y, C_vec.z, res_i_occ, res_i_b_fac );
			res_i->newAtom(CB, CB_vec.x, CB_vec.y, CB_vec.z, res_i_occ, res_i_b_fac);
			res_i->newAtom(O, O_vec.x, O_vec.y, O_vec.z, res_i_occ, res_i_b_fac);
			res_iplus1->newAtom(N, N_vec.x, N_vec.y, N_vec.z, res_iplus1_occ, res_iplus1_b_fac);
			 */
		}
		else {
			int used_key;
			if (withNoise){
				used_key = this->PRO_findClosestKey_with_noise(overall_key, noise_level);
			}
			else {
				used_key = PRO_findClosestKey(overall_key);
			}

			vector3d C_local;
			vector3d CB_local;
			vector3d O_local;
			vector3d N_local;

			PRO_data_map.find(used_key)->second.get_local_vecs(
					C_local,
					O_local,
					N_local,
					CB_local,
					first_dih_bin,
					last_dih_bin);

			const vector3d C_vec = ele.CA_atoms[1]->get_coords() + (C_local.x * x_unit
					+ C_local.y * y_unit
					+ C_local.z * z_unit);
			const vector3d CB_vec = ele.CA_atoms[1]->get_coords() + (CB_local.x * x_unit
					+ CB_local.y * y_unit
					+ CB_local.z * z_unit);
			const vector3d O_vec = ele.CA_atoms[1]->get_coords() + (O_local.x * x_unit
					+ O_local.y * y_unit
					+ O_local.z * z_unit);
			const vector3d N_vec = ele.CA_atoms[1]->get_coords() + (N_local.x * x_unit
					+ N_local.y * y_unit
					+ N_local.z * z_unit);

			if (!pose_->get_bb_atom(POSE::C, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "C" << endl;
			}
			pose_->add_new_atom(C_vec, ref_C, resNum_1);
			if (!(pose_->get_residue(resNum_1)->get_type() == residue_type("GLY")) && changeCB){
				if (!pose_->get_bb_atom(POSE::CB, resNum_1)->isSet()){
					cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
							<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
							<< "CB" << endl;
				}
				pose_->add_new_atom(CB_vec, ref_CB, resNum_1);
			}
			if (!pose_->get_bb_atom(POSE::O, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "O" << endl;
			}
			pose_->add_new_atom(O_vec, ref_O, resNum_1);
			if (!pose_->get_bb_atom(POSE::N, resNum_1+1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1+1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1+1)->get_chain()->getChainID() << "\t"
						<< "N" << endl;
			}
			pose_->add_new_atom(N_vec, ref_N, resNum_1+1);

			/*
			res_i->newAtom(C, C_vec.x, C_vec.y, C_vec.z, res_i_occ, res_i_b_fac);
			res_i->newAtom(CB, CB_vec.x, CB_vec.y, CB_vec.z, res_i_occ, res_i_b_fac);
			res_i->newAtom(O, O_vec.x, O_vec.y, O_vec.z, res_i_occ, res_i_b_fac);
			res_iplus1->newAtom(N, N_vec.x, N_vec.y, N_vec.z, res_iplus1_occ, res_iplus1_b_fac);
			 */
		}
		/*
		res_i->resetCBcoor();
		res_iplus1->resetCBcoor();
		res_iplus1->addH();
		 */


	}
	else {
		// IBBB fall back


		const int overall_key = ele.IBBB_key;


		//const vector3d v_i_minus1 = ele.CA_atoms[1]->get_coords() - ele.CA_atoms[0]->get_coords();
		const vector3d v_i = ele.CA_atoms[2]->get_coords() - ele.CA_atoms[1]->get_coords();
		//const vector3d v_i_plus1 = ele.CA_atoms[3]->get_coords() - ele.CA_atoms[2]->get_coords();
		const vector3d v_p = ele.CA_atoms[3]->get_coords() - ele.CA_atoms[1]->get_coords();

		const vector3d x_unit = (v_i * v_p) / (v_i * v_p).mod();
		const vector3d y_unit = (v_p * x_unit) / (v_p * x_unit).mod();
		const vector3d z_unit = ( x_unit * y_unit ) / ( x_unit * y_unit ).mod();

		if (!(pose_->get_residue(resNum_1)->get_type() == residue_type("PRO"))){

			int used_key;
			if (withNoise){
				used_key = this->IBBBfindClosestKey_with_noise(overall_key, noise_level);
			}
			else {
				used_key = IBBBfindClosestKey(overall_key);
			}

			vector3d C_local = C_vec_map.find(used_key) !=  C_vec_map.end() ? C_vec_map.find(used_key)->second : vector3d();
			vector3d CB_local = CB_vec_map.find(used_key) !=  CB_vec_map.end() ? CB_vec_map.find(used_key)->second : vector3d();
			vector3d O_local = O_vec_map.find(used_key) !=  O_vec_map.end() ? O_vec_map.find(used_key)->second : vector3d();
			vector3d N_local = N_vec_map.find(used_key) !=  N_vec_map.end() ? N_vec_map.find(used_key)->second : vector3d();

			/*
			non_PRO_data_map.find(used_key)->second.get_local_vecs(
					C_local,
					O_local,
					N_local,
					CB_local,
					first_dih_bin,
					last_dih_bin);
					*/

			const vector3d C_vec = ele.CA_atoms[1]->get_coords() + (C_local.x * x_unit
					+ C_local.y * y_unit
					+ C_local.z * z_unit);
			const vector3d CB_vec = ele.CA_atoms[1]->get_coords() + (CB_local.x * x_unit
					+ CB_local.y * y_unit
					+ CB_local.z * z_unit);
			const vector3d O_vec = ele.CA_atoms[1]->get_coords() + (O_local.x * x_unit
					+ O_local.y * y_unit
					+ O_local.z * z_unit);
			const vector3d N_vec = ele.CA_atoms[1]->get_coords() + (N_local.x * x_unit
					+ N_local.y * y_unit
					+ N_local.z * z_unit);


			if (!pose_->get_bb_atom(POSE::C, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "C" << endl;
			}
			pose_->add_new_atom(C_vec, ref_C, resNum_1);
			if (!(pose_->get_residue(resNum_1)->get_type() == residue_type("GLY")) && changeCB){
				if (!pose_->get_bb_atom(POSE::CB, resNum_1)->isSet()){
					cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
							<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
							<< "CB" << endl;
				}
				pose_->add_new_atom(CB_vec, ref_CB, resNum_1);
			}
			if (!pose_->get_bb_atom(POSE::O, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "O" << endl;
			}
			pose_->add_new_atom(O_vec, ref_O, resNum_1);
			if (!pose_->get_bb_atom(POSE::N, resNum_1+1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1+1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1+1)->get_chain()->getChainID() << "\t"
						<< "N" << endl;
			}
			pose_->add_new_atom(N_vec, ref_N, resNum_1+1);
			/*
			pose_->add_new_atom(C_vec, ref_C, resNum_1);
			pose_->add_new_atom(CB_vec, ref_CB, resNum_1);
			pose_->add_new_atom(O_vec, ref_O, resNum_1);
			pose_->add_new_atom(N_vec, ref_N, resNum_1+1);
			 */

			/*
			res_i->newAtom(C, C_vec.x, C_vec.y, C_vec.z, res_i_occ, res_i_b_fac );
			res_i->newAtom(CB, CB_vec.x, CB_vec.y, CB_vec.z, res_i_occ, res_i_b_fac);
			res_i->newAtom(O, O_vec.x, O_vec.y, O_vec.z, res_i_occ, res_i_b_fac);
			res_iplus1->newAtom(N, N_vec.x, N_vec.y, N_vec.z, res_iplus1_occ, res_iplus1_b_fac);
			 */
		}
		else {
			int used_key;
			if (withNoise){
				used_key = this->PRO_IBBBfindClosestKey_with_noise(overall_key, noise_level);
			}
			else {
				used_key = PRO_IBBBfindClosestKey(overall_key);
			}


			vector3d C_local = PRO_C_vec_map.find(used_key) !=  PRO_C_vec_map.end() ? PRO_C_vec_map.find(used_key)->second : vector3d();
			vector3d CB_local = PRO_CB_vec_map.find(used_key) !=  PRO_CB_vec_map.end() ? PRO_CB_vec_map.find(used_key)->second : vector3d();
			vector3d O_local = PRO_O_vec_map.find(used_key) !=  PRO_O_vec_map.end() ? PRO_O_vec_map.find(used_key)->second : vector3d();
			vector3d N_local = PRO_N_vec_map.find(used_key) !=  PRO_N_vec_map.end() ? PRO_N_vec_map.find(used_key)->second : vector3d();

			/*
			PRO_data_map.find(used_key)->second.get_local_vecs(
					C_local,
					O_local,
					N_local,
					CB_local,
					first_dih_bin,
					last_dih_bin);
			 */

			const vector3d C_vec = ele.CA_atoms[1]->get_coords() + (C_local.x * x_unit
					+ C_local.y * y_unit
					+ C_local.z * z_unit);
			const vector3d CB_vec = ele.CA_atoms[1]->get_coords() + (CB_local.x * x_unit
					+ CB_local.y * y_unit
					+ CB_local.z * z_unit);
			const vector3d O_vec = ele.CA_atoms[1]->get_coords() + (O_local.x * x_unit
					+ O_local.y * y_unit
					+ O_local.z * z_unit);
			const vector3d N_vec = ele.CA_atoms[1]->get_coords() + (N_local.x * x_unit
					+ N_local.y * y_unit
					+ N_local.z * z_unit);

			if (!pose_->get_bb_atom(POSE::C, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "C" << endl;
			}
			pose_->add_new_atom(C_vec, ref_C, resNum_1);
			if (!(pose_->get_residue(resNum_1)->get_type() == residue_type("GLY")) && changeCB){
				if (!pose_->get_bb_atom(POSE::CB, resNum_1)->isSet()){
					cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
							<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
							<< "CB" << endl;
				}
				pose_->add_new_atom(CB_vec, ref_CB, resNum_1);
			}
			if (!pose_->get_bb_atom(POSE::O, resNum_1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1)->get_chain()->getChainID() << "\t"
						<< "O" << endl;
			}
			pose_->add_new_atom(O_vec, ref_O, resNum_1);
			if (!pose_->get_bb_atom(POSE::N, resNum_1+1)->isSet()){
				cout << "backbone_builder: adding new atom at: " << pose_->get_residue(resNum_1+1)->get_pdb_residue_index() << "\t"
						<< pose_->get_residue(resNum_1+1)->get_chain()->getChainID() << "\t"
						<< "N" << endl;
			}
			pose_->add_new_atom(N_vec, ref_N, resNum_1+1);
			/*
			pose_->add_new_atom(C_vec, ref_C, resNum_1);
			pose_->add_new_atom(CB_vec, ref_CB, resNum_1);
			pose_->add_new_atom(O_vec, ref_O, resNum_1);
			pose_->add_new_atom(N_vec, ref_N, resNum_1+1);
			*/

		}


	}

	// TODO use better method
	POSE_UTILS::quick_add_HN(pose_,resNum_1+1, true);

	return true;

}

int backbone_builder::findClosestKey(const int search_key) const{
	int_bb_ele_umap::const_iterator iter;
	double best_dist = FLT_MAX;
	int best_key = 0;
	if (non_PRO_data_map.find(search_key) != non_PRO_data_map.end()){
		return search_key;
	}
	for ( iter = non_PRO_data_map.begin(); iter != non_PRO_data_map.end(); iter++ ){
		int key = iter->first;
		double dist = key_distance(search_key, key);

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}
	return best_key;
}

int backbone_builder::PRO_findClosestKey(const int search_key) const{
	int_bb_ele_umap::const_iterator iter;
	double best_dist = FLT_MAX;
	int best_key = 0;
	if (PRO_data_map.find(search_key) != PRO_data_map.end()){
		return search_key;
	}
	for ( iter = PRO_data_map.begin(); iter != PRO_data_map.end(); iter++ ){
		int key = iter->first;
		double dist = key_distance(search_key, key);

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}
	return best_key;
}

int backbone_builder::findClosestKey_with_noise(const int search_key, const double noise_dist) const{
	int_bb_ele_umap::const_iterator iter;
	int_vector candidates;
	double best_dist = FLT_MAX;
	int best_key = 0;


	for ( iter = non_PRO_data_map.begin(); iter != non_PRO_data_map.end(); iter++ ){
		int key = iter->first;
		double dist = key_distance(search_key, key);

		if (dist < noise_dist ){
			candidates.push_back(key);
		}

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}

	const int list_size = static_cast<int>(candidates.size());
	if (list_size > 1){
		best_key = candidates[rand_gen->randInt(list_size - 1)];
	}



	return best_key;
}

int backbone_builder::PRO_findClosestKey_with_noise(const int search_key, const double noise_dist) const{
	int_bb_ele_umap::const_iterator iter;
	int_vector candidates;
	double best_dist = FLT_MAX;
	int best_key = 0;


	for ( iter = PRO_data_map.begin(); iter != PRO_data_map.end(); iter++ ){
		int key = iter->first;
		double dist = key_distance(search_key, key);

		if (dist < noise_dist ){
			candidates.push_back(key);
		}

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}

	const int list_size = static_cast<int>(candidates.size());
	if (list_size > 1){
		best_key = candidates[rand_gen->randInt(list_size - 1)];
	}



	return best_key;
}





ostream &operator<<( ostream& output,  backbone_builder & db ){

	int_bb_ele_umap::const_iterator iter;
	int_bb_dih_ele_map::const_iterator dih_iter;

	output << "DIST_BIN_SIZE\t" << db.bin_size << endl;

	for ( iter = db.non_PRO_data_map.begin(); iter != db.non_PRO_data_map.end(); iter++ ){
		for (dih_iter = iter->second.dih_ele_map.begin(); dih_iter != iter->second.dih_ele_map.end(); dih_iter++){

			output << iter->first << "\t"
					<< dih_iter->first << "\t"
					<< dih_iter->second.count << "\t"
					<< "REG" << "\t"
					<< (dih_iter->second.C_local)
					<< (dih_iter->second.O_local)
					<< (dih_iter->second.N_local)
					<< (dih_iter->second.CB_local)
					<< endl;
		}
	}

	for ( iter = db.PRO_data_map.begin(); iter != db.PRO_data_map.end(); iter++ ){
		for (dih_iter = iter->second.dih_ele_map.begin(); dih_iter != iter->second.dih_ele_map.end(); dih_iter++){

			output << iter->first << "\t"
					<< dih_iter->first << "\t"
					<< dih_iter->second.count << "\t"
					<< "PRO" << "\t"
					<< (dih_iter->second.C_local)
					<< (dih_iter->second.O_local)
					<< (dih_iter->second.N_local)
					<< (dih_iter->second.CB_local)
					<< endl;
		}
	}

	return output;
}

std::istream& backbone_builder::load_bb_db( std::istream& input){
	string lineStr;
	long length, lineNum = 0 ;


	non_PRO_data_map.clear();
	PRO_data_map.clear();


	split_vector_type SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			const string paraName = SplitVec[0];
			if (paraName.compare("DIST_BIN_SIZE") != 0){
				const int main_key = lexical_cast<int>( SplitVec[0] );
				const int dih_key = lexical_cast<int>( SplitVec[1] );
				long freq = lexical_cast<long>( SplitVec[2] );

				if (SplitVec[3].compare("REG") == 0 ) {

					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].count = static_cast<double>(freq);

					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.x = lexical_cast<double>(SplitVec[4]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.y = lexical_cast<double>(SplitVec[5]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.z = lexical_cast<double>(SplitVec[6]);

					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.x = lexical_cast<double>(SplitVec[7]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.y = lexical_cast<double>(SplitVec[8]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.z = lexical_cast<double>(SplitVec[9]);

					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.x = lexical_cast<double>(SplitVec[10]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.y = lexical_cast<double>(SplitVec[11]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.z = lexical_cast<double>(SplitVec[12]);

					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.x = lexical_cast<double>(SplitVec[13]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.y = lexical_cast<double>(SplitVec[14]);
					(non_PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.z = lexical_cast<double>(SplitVec[15]);
				}
				else if (SplitVec[3].compare("PRO") == 0 ) {
					(PRO_data_map[main_key]).dih_ele_map[dih_key].count = static_cast<double>(freq);

					(PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.x = lexical_cast<double>(SplitVec[4]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.y = lexical_cast<double>(SplitVec[5]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.z = lexical_cast<double>(SplitVec[6]);

					(PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.x = lexical_cast<double>(SplitVec[7]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.y = lexical_cast<double>(SplitVec[8]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.z = lexical_cast<double>(SplitVec[9]);

					(PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.x = lexical_cast<double>(SplitVec[10]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.y = lexical_cast<double>(SplitVec[11]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.z = lexical_cast<double>(SplitVec[12]);

					(PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.x = lexical_cast<double>(SplitVec[13]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.y = lexical_cast<double>(SplitVec[14]);
					(PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.z = lexical_cast<double>(SplitVec[15]);
				}

			}
			else {
				const double this_bin_size = lexical_cast<double>( SplitVec[1] );

				bin_size = this_bin_size;

			}
		}
	}


	return input;
}

std::istream& backbone_builder::load_bb_IBBB_db( std::istream& input){
	string lineStr;
	long length, lineNum = 0 ;


	int_overall_bin_count.clear();
	C_vec_map.clear();
	O_vec_map.clear();
	N_vec_map.clear();

	PRO_int_overall_bin_count.clear();
	PRO_C_vec_map.clear();
	PRO_O_vec_map.clear();
	PRO_N_vec_map.clear();

	split_vector_type SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			int key = lexical_cast<int>( SplitVec[0] );
			long freq = lexical_cast<long>( SplitVec[1] );
			if (SplitVec[2].compare("REG") == 0 ) {

				(int_overall_bin_count[key]) = freq;

				(C_vec_map[key]).x = lexical_cast<double>(SplitVec[3]);
				(C_vec_map[key]).y = lexical_cast<double>(SplitVec[4]);
				(C_vec_map[key]).z = lexical_cast<double>(SplitVec[5]);

				(O_vec_map[key]).x = lexical_cast<double>(SplitVec[6]);
				(O_vec_map[key]).y = lexical_cast<double>(SplitVec[7]);
				(O_vec_map[key]).z = lexical_cast<double>(SplitVec[8]);

				(N_vec_map[key]).x = lexical_cast<double>(SplitVec[9]);
				(N_vec_map[key]).y = lexical_cast<double>(SplitVec[10]);
				(N_vec_map[key]).z = lexical_cast<double>(SplitVec[11]);

				(CB_vec_map[key]).x = lexical_cast<double>(SplitVec[12]);
				(CB_vec_map[key]).y = lexical_cast<double>(SplitVec[13]);
				(CB_vec_map[key]).z = lexical_cast<double>(SplitVec[14]);
			}
			else if (SplitVec[2].compare("PRO") == 0 ) {
				(PRO_int_overall_bin_count[key]) = freq;

				(PRO_C_vec_map[key]).x = lexical_cast<double>(SplitVec[3]);
				(PRO_C_vec_map[key]).y = lexical_cast<double>(SplitVec[4]);
				(PRO_C_vec_map[key]).z = lexical_cast<double>(SplitVec[5]);

				(PRO_O_vec_map[key]).x = lexical_cast<double>(SplitVec[6]);
				(PRO_O_vec_map[key]).y = lexical_cast<double>(SplitVec[7]);
				(PRO_O_vec_map[key]).z = lexical_cast<double>(SplitVec[8]);

				(PRO_N_vec_map[key]).x = lexical_cast<double>(SplitVec[9]);
				(PRO_N_vec_map[key]).y = lexical_cast<double>(SplitVec[10]);
				(PRO_N_vec_map[key]).z = lexical_cast<double>(SplitVec[11]);

				(PRO_CB_vec_map[key]).x = lexical_cast<double>(SplitVec[12]);
				(PRO_CB_vec_map[key]).y = lexical_cast<double>(SplitVec[13]);
				(PRO_CB_vec_map[key]).z = lexical_cast<double>(SplitVec[14]);
			}

		}
	}


	return input;
}

/*
istream &operator>>( istream& input, backbone_builder & db ) {


	string lineStr;
	long length, lineNum = 0 ;


	db.non_PRO_data_map.clear();
	db.PRO_data_map.clear();


	split_vector_type SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t") );
			const string paraName = SplitVec[0];
			if (paraName.compare("DIST_BIN_SIZE") != 0){
				const int main_key = lexical_cast<int>( SplitVec[0] );
				const int dih_key = lexical_cast<int>( SplitVec[1] );
				long freq = lexical_cast<long>( SplitVec[2] );

				if (SplitVec[3].compare("REG") == 0 ) {

					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].count = static_cast<double>(freq);

					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.x = lexical_cast<double>(SplitVec[4]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.y = lexical_cast<double>(SplitVec[5]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.z = lexical_cast<double>(SplitVec[6]);

					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.x = lexical_cast<double>(SplitVec[7]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.y = lexical_cast<double>(SplitVec[8]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.z = lexical_cast<double>(SplitVec[9]);

					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.x = lexical_cast<double>(SplitVec[10]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.y = lexical_cast<double>(SplitVec[11]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.z = lexical_cast<double>(SplitVec[12]);

					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.x = lexical_cast<double>(SplitVec[13]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.y = lexical_cast<double>(SplitVec[14]);
					(db.non_PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.z = lexical_cast<double>(SplitVec[15]);
				}
				else if (SplitVec[3].compare("PRO") == 0 ) {
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].count = static_cast<double>(freq);

					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.x = lexical_cast<double>(SplitVec[4]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.y = lexical_cast<double>(SplitVec[5]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].C_local.z = lexical_cast<double>(SplitVec[6]);

					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.x = lexical_cast<double>(SplitVec[7]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.y = lexical_cast<double>(SplitVec[8]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].O_local.z = lexical_cast<double>(SplitVec[9]);

					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.x = lexical_cast<double>(SplitVec[10]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.y = lexical_cast<double>(SplitVec[11]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].N_local.z = lexical_cast<double>(SplitVec[12]);

					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.x = lexical_cast<double>(SplitVec[13]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.y = lexical_cast<double>(SplitVec[14]);
					(db.PRO_data_map[main_key]).dih_ele_map[dih_key].CB_local.z = lexical_cast<double>(SplitVec[15]);
				}

			}
			else {
				const double bin_size = lexical_cast<double>( SplitVec[1] );

				db.bin_size = bin_size;

			}
		}
	}


	return input;
}
*/


double backbone_builder::key_distance(const int key1, const int key2) const{

	//const int abs_key1;// = abs(key1);
	//const int abs_key2;// = abs(key2);

	int sign_key1;// = key1 / abs_key1;
	int sign_key2 ;//= key2 / abs_key2;

	int key1_d_m1_p2_bin;// = sign_key1 * (abs_key1 / 10000);
	int key1_d_0_p2_bin;// = (abs_key1 - abs(key1_d_m1_p2_bin * 10000))
						//		/ 100;
	int key1_d_m1_p1_bin;// = (abs_key1 - abs(key1_d_m1_p2_bin * 10000)
						 //		- (key1_d_0_p2_bin * 100));

	int key1_first_dih_bin;
	int key1_last_dih_bin;

	this->decomposeKey(key1,
			sign_key1,
			key1_d_m1_p2_bin,
			key1_d_0_p2_bin,
			key1_d_m1_p1_bin,
			key1_first_dih_bin,
			key1_last_dih_bin);

	key1_d_m1_p2_bin = key1_d_m1_p2_bin * sign_key1;


	int key2_d_m1_p2_bin;// = sign_key2 * (abs_key2 / 10000);
	int key2_d_0_p2_bin;// = (abs_key2 - abs(key2_d_m1_p2_bin * 10000))
						//		/ 100;
	int key2_d_m1_p1_bin;// = (abs_key2 - abs(key2_d_m1_p2_bin * 10000)
						 //		- (key2_d_0_p2_bin * 100));

	int key2_first_dih_bin;
	int key2_last_dih_bin;

	this->decomposeKey(key2,
			sign_key2,
			key2_d_m1_p2_bin,
			key2_d_0_p2_bin,
			key2_d_m1_p1_bin,
			key2_first_dih_bin,
			key2_last_dih_bin);
	key2_d_m1_p2_bin = key2_d_m1_p2_bin * sign_key2;

	/*
	cout << key1 << "\t"
		 << key1_d_m1_p2_bin << "\t"
		 << key1_d_0_p2_bin << "\t"
		 << key1_d_m1_p1_bin << "\t"
		 << endl;

	cout << key2 << "\t"
		 << key2_d_m1_p2_bin << "\t"
		 << key2_d_0_p2_bin << "\t"
		 << key2_d_m1_p1_bin << "\t"
		 << endl;
	*/
	double distance = sqrt(static_cast<double>((key1_d_m1_p2_bin - key2_d_m1_p2_bin) * (key1_d_m1_p2_bin - key2_d_m1_p2_bin))
					+ static_cast<double>((key1_d_0_p2_bin - key2_d_0_p2_bin) * (key1_d_0_p2_bin - key2_d_0_p2_bin))
					+ static_cast<double>((key1_d_m1_p1_bin - key2_d_m1_p1_bin) * (key1_d_m1_p1_bin - key2_d_m1_p1_bin))
					+ dih_dist_weight * static_cast<double>((key1_first_dih_bin - key2_first_dih_bin) * (key1_first_dih_bin - key2_first_dih_bin))
					+ dih_dist_weight * static_cast<double>((key1_last_dih_bin - key2_last_dih_bin) * (key1_last_dih_bin - key2_last_dih_bin)));

	/*
	cout << distance << endl;
	*/

	return distance;
}

double backbone_builder::IBBB_key_distance(const int key1, const int key2) const{

	const int abs_key1 = abs(key1);
	const int abs_key2 = abs(key2);

	const int sign_key1 = key1 / abs_key1;
	const int sign_key2 = key2 / abs_key2;

	/*
	const int overall_key = static_cast<int>(chi) * (10000 * abs(d_m1_p2_bin)
							+ 100 * d_0_p2_bin
							+  d_m1_p1_bin);
	*/
	const int key1_d_m1_p2_bin = sign_key1 * (abs_key1 / 10000);
	const int key1_d_0_p2_bin = (abs_key1 - abs(key1_d_m1_p2_bin * 10000))
								/ 100;
	const int key1_d_m1_p1_bin = (abs_key1 - abs(key1_d_m1_p2_bin * 10000)
								- (key1_d_0_p2_bin * 100));

	const int key2_d_m1_p2_bin = sign_key2 * (abs_key2 / 10000);
	const int key2_d_0_p2_bin = (abs_key2 - abs(key2_d_m1_p2_bin * 10000))
								/ 100;
	const int key2_d_m1_p1_bin = (abs_key2 - abs(key2_d_m1_p2_bin * 10000)
								- (key2_d_0_p2_bin * 100));

	/*
	cout << key1 << "\t"
		 << key1_d_m1_p2_bin << "\t"
		 << key1_d_0_p2_bin << "\t"
		 << key1_d_m1_p1_bin << "\t"
		 << endl;

	cout << key2 << "\t"
		 << key2_d_m1_p2_bin << "\t"
		 << key2_d_0_p2_bin << "\t"
		 << key2_d_m1_p1_bin << "\t"
		 << endl;
	*/
	double distance = sqrt((key1_d_m1_p2_bin - key2_d_m1_p2_bin) * (key1_d_m1_p2_bin - key2_d_m1_p2_bin)
					+ (key1_d_0_p2_bin - key2_d_0_p2_bin) * (key1_d_0_p2_bin - key2_d_0_p2_bin)
					+ (key1_d_m1_p1_bin - key2_d_m1_p1_bin) * (key1_d_m1_p1_bin - key2_d_m1_p1_bin));

	/*
	cout << distance << endl;
	*/

	return distance;
}

// TODO implement this
bool backbone_builder::validate_residues(
		PRODART::POSE::const_residue_shared_ptr res_i_minus2,
		PRODART::POSE::const_residue_shared_ptr res_i_minus1, PRODART::POSE::const_residue_shared_ptr res_i,
		PRODART::POSE::const_residue_shared_ptr res_iplus1, PRODART::POSE::const_residue_shared_ptr res_iplus2,
		PRODART::POSE::const_residue_shared_ptr res_iplus3
		){
	/*
	if (res_i_minus2 != 0 && res_i_minus1 != 0 && res_i != 0
			&& res_iplus1 != 0 && res_iplus2 != 0 &&
			res_iplus3 != 0){

		if (	res_i_minus2->getChain() == res_i_minus1->getChain() &&
				res_i_minus1->getChain() == res_i->getChain() &&
				res_i->getChain() == res_iplus1->getChain() &&
				res_iplus1->getChain() == res_iplus2->getChain() &&
				res_iplus2->getChain() == res_iplus3->getChain()
				) {

			if (	res_i_minus2->hasFullBB()
					&& res_i_minus1->hasFullBB() && res_i->hasFullBB()
					&& res_iplus1->hasFullBB() && res_iplus2->hasFullBB()
					&& res_iplus3->hasFullBB()
					){

				vector3d CA_m2 = res_i_minus2->getAtom(CA)->getVec();
				vector3d CA_m1 = res_i_minus1->getAtom(CA)->getVec();
				vector3d CA_0 = res_i->getAtom(CA)->getVec();
				vector3d CA_p1 = res_iplus1->getAtom(CA)->getVec();
				vector3d CA_p2 = res_iplus2->getAtom(CA)->getVec();
				vector3d CA_p3 = res_iplus3->getAtom(CA)->getVec();

				if ( 	fabs((CA_m2 - CA_m1).mod() - 3.4) < 0.8 &&
						fabs((CA_m1 - CA_0).mod() - 3.4) < 0.8 &&
						fabs((CA_0 - CA_p1).mod() - 3.4) < 0.8 &&
						fabs((CA_p1 - CA_p2).mod() - 3.4) < 0.8 &&
						fabs((CA_p2 - CA_p3).mod() - 3.4) < 0.8){

					vector3d C_0 = res_i->getAtom(C)->getVec();
					vector3d N_p1 = res_iplus1->getAtom(N)->getVec();

					if (fabs((CA_0 - C_0).mod() - 1.51) < 0.07 &&
						fabs((C_0 - N_p1).mod() - 1.33) < 0.07 &&
						fabs((N_p1 - CA_p1).mod() - 1.46) < 0.07 ){

						return true;
					}
				}
			}


		}


	}
	//cerr << "ERROR: IBBB_validate_residues: residues not validated\n";
	return false;

	*/
	return true;
}





int backbone_builder::IBBBfindClosestKey(const int search_key) const{
	int_ulong_umap::const_iterator iter;
	double best_dist = FLT_MAX;
	int best_key = 0;
	if (int_overall_bin_count.find(search_key) != int_overall_bin_count.end()){
		return search_key;
	}
	for ( iter = int_overall_bin_count.begin(); iter != int_overall_bin_count.end(); iter++ ){
		int key = iter->first;
		double dist = IBBB_key_distance(search_key, key);

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}
	return best_key;
}

int backbone_builder::PRO_IBBBfindClosestKey(const int search_key) const{
	int_ulong_umap::const_iterator iter;
	double best_dist = FLT_MAX;
	int best_key = 0;
	if (PRO_int_overall_bin_count.find(search_key) != PRO_int_overall_bin_count.end()){
		return search_key;
	}
	for ( iter = PRO_int_overall_bin_count.begin(); iter != PRO_int_overall_bin_count.end(); iter++ ){
		int key = iter->first;
		double dist = IBBB_key_distance(search_key, key);

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}
	return best_key;
}


int backbone_builder::IBBBfindClosestKey_with_noise(const int search_key, const double noise_dist) const{
	int_ulong_umap::const_iterator iter;
	int_vector candidates;
	double best_dist = FLT_MAX;
	int best_key = 0;


	for ( iter = int_overall_bin_count.begin(); iter != int_overall_bin_count.end(); iter++ ){
		int key = iter->first;
		double dist = IBBB_key_distance(search_key, key);

		if (dist < noise_dist ){
			candidates.push_back(key);
		}

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}

	const int list_size = static_cast<int>(candidates.size());
	if (list_size > 1){
		best_key = candidates[rand_gen->randInt(list_size - 1)];
	}



	return best_key;
}

int backbone_builder::PRO_IBBBfindClosestKey_with_noise(const int search_key, const double noise_dist) const{
	int_ulong_umap::const_iterator iter;
	int_vector candidates;
	double best_dist = FLT_MAX;
	int best_key = 0;


	for ( iter = PRO_int_overall_bin_count.begin(); iter != PRO_int_overall_bin_count.end(); iter++ ){
		int key = iter->first;
		double dist = IBBB_key_distance(search_key, key);

		if (dist < noise_dist ){
			candidates.push_back(key);
		}

		if (dist < best_dist ){
			best_dist = dist;
			best_key = key;
		}

	}

	const int list_size = static_cast<int>(candidates.size());
	if (list_size > 1){
		best_key = candidates[rand_gen->randInt(list_size - 1)];
	}



	return best_key;
}


}
}
}

