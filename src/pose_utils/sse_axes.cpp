/*
 * sse_axes.cpp
 *
 *  Created on: Mar 12, 2012
 *      Author: jmacdon
 */
#include "sse_axes.h"

using namespace PRODART::POSE;
using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::pow;
using std::sqrt;
using std::asin;
using std::string;

using boost::lexical_cast;
using boost::split;
using boost::is_any_of;
using boost::trim;

using std::ios;
using std::setiosflags;
using std::resetiosflags;
using std::setprecision;
using std::setw;

using PRODART::UTILS::vector3d;
using namespace PRODART::UTILS;

namespace PRODART {
namespace POSE_UTILS {


void calculate_helix_axis(const PRODART::UTILS::vector3d& c1,
		const PRODART::UTILS::vector3d& c2,
		const PRODART::UTILS::vector3d& c3,
		const PRODART::UTILS::vector3d& c4,
		PRODART::UTILS::vector3d& start,
		PRODART::UTILS::vector3d& end){
	PRODART::UTILS::vector3d v1 = (c1 - c2) + (c3 - c2);
	v1 = v1 / v1.mod();

	PRODART::UTILS::vector3d v2 = (c2 - c3) + (c4 - c3);
	v2 = v2 / v2.mod();

	PRODART::UTILS::vector3d H = v1 * v2;
	H = H / H.mod();

	const double d = fabs((c3 - c2).dot(H));

	const double radius = fabs(( (H * d).mod_sq() - (c3 - c2).mod_sq()) /
			(2.0 * ( fabs((c2 - c3).dot(v2)) )));

	//std::cout << "radius:\t" << radius << endl;

	start = c2 + (radius * v1);
	end = c3 + (radius * v2);
}

void calculate_helix_axis_long(const vector3d& c1,
		const vector3d& c2,
		const vector3d& c3,
		const vector3d& c4,
		vector3d& start,
		vector3d& end){
	vector3d v1 = (c1 - c2) + (c3 - c2);
	v1 = v1 / v1.mod();

	vector3d v2 = (c2 - c3) + (c4 - c3);
	v2 = v2 / v2.mod();

	vector3d H = v1 * v2;
	H = H / H.mod();

	const double d = fabs((c3 - c2).dot(H));

	const double radius = fabs(( (H * d).mod_sq() - (c3 - c2).mod_sq()) /
			(2.0 * ( fabs((c2 - c3).dot(v2)) )));

	//std::cout << "radius:\t" << radius << endl;

	const vector3d short_start = c2 + (radius * v1);
	const vector3d short_end = c3 + (radius * v2);

	const vector3d axis_vec = (short_end - short_start) / (short_end - short_start).mod();


	const vector3d first = c1 - short_start;
	const vector3d last = c4 - short_start;

	start = short_start + (first.dot(axis_vec) * axis_vec);
	end = short_start + (last.dot(axis_vec) * axis_vec);


}

void calculate_strand_axis(const PRODART::UTILS::vector3d& c1,
		const PRODART::UTILS::vector3d& c2,
		const PRODART::UTILS::vector3d& c3,
		const PRODART::UTILS::vector3d& c4,
		PRODART::UTILS::vector3d& start,
		PRODART::UTILS::vector3d& end){

	start = (c1 + c2 + c3) / 3.0;
	end = (c2 + c3 + c4) / 3.0;

}

void calculate_strand_axis_long(const vector3d& c1,
		const vector3d& c2,
		const vector3d& c3,
		const vector3d& c4,
		vector3d& start,
		vector3d& end){

	const vector3d short_start = (c1 + c2 + c3) / 3.0;
	const vector3d short_end = (c2 + c3 + c4) / 3.0;

	const vector3d axis_vec = (short_end - short_start) / (short_end - short_start).mod();


	const vector3d first = c1 - short_start;
	const vector3d last = c4 - short_start;

	start = short_start + (first.dot(axis_vec) * axis_vec);
	end = short_start + (last.dot(axis_vec) * axis_vec);

}

bool get_helix_end_points_long(const PRODART::POSE::const_pose_shared_ptr protein,
		const int start_index,
		const int end_index,
		vector3d& startVec,
		vector3d& endVec){

	if (end_index - start_index <= 2){
		return false;
	}

	PRODART::UTILS::vector3d_vector coords;


	for (int i = start_index; i <= end_index; i++){
		coords.push_back(protein->get_bb_atom_coords(POSE::CA,i));
	}

	return get_helix_end_points_long(coords, startVec, endVec);

	/*
	else if (end_index - start_index <= 3){

		//protein->get_bb_atom_coords(POSE::CA, 0);
		vector3d c1 = protein->get_bb_atom_coords(POSE::CA,start_index);
		vector3d c2 = protein->get_bb_atom_coords(POSE::CA,start_index+1);
		vector3d c3 = protein->get_bb_atom_coords(POSE::CA,start_index+2);
		vector3d c4 = protein->get_bb_atom_coords(POSE::CA,start_index+3);

		calculate_helix_axis_long(c1,c2,c3,c4,startVec,endVec);

		return true;

	}

	PRODART::UTILS::vector3d_vector coords;

	coords.reserve(end_index - start_index);

	vector3d c1 = protein->get_bb_atom_coords(POSE::CA,start_index);
	vector3d c2 = protein->get_bb_atom_coords(POSE::CA,start_index+1);
	vector3d c3 = protein->get_bb_atom_coords(POSE::CA,start_index+2);
	vector3d c4 = protein->get_bb_atom_coords(POSE::CA,start_index+3);

	for (int i = start_index+4; i <= end_index; i++){

		vector3d  a,b;

		calculate_helix_axis(c1,c2,c3,c4,a,b);
		coords.push_back(a);
		coords.push_back(b);

		c1 = c2;
		c2 = c3;
		c3 = c4;
		c4 = protein->get_bb_atom_coords(POSE::CA,i);


	}

	const vector3d first_atom = protein->get_bb_atom_coords(POSE::CA,start_index);
	const vector3d last_atom = protein->get_bb_atom_coords(POSE::CA,end_index);
	PRODART::UTILS::line_fit3d(coords,
			startVec,
			endVec,
			first_atom,
			last_atom);

	return true;
	*/

}


bool get_strand_end_points_long(const PRODART::POSE::const_pose_shared_ptr protein,
		const int start_index,
		const int end_index,
		vector3d& startVec,
		vector3d& endVec){

	if (end_index - start_index <= 2){
		return false;
	}

	PRODART::UTILS::vector3d_vector coords;


	for (int i = start_index; i <= end_index; i++){
		coords.push_back(protein->get_bb_atom_coords(POSE::CA,i));
	}

	return get_strand_end_points_long(coords, startVec, endVec);

	/*
	else if (end_index - start_index <= 3){

		vector3d c1 = protein->get_bb_atom_coords(POSE::CA,start_index);
		vector3d c2 = protein->get_bb_atom_coords(POSE::CA,start_index+1);
		vector3d c3 = protein->get_bb_atom_coords(POSE::CA,start_index+2);
		vector3d c4 = protein->get_bb_atom_coords(POSE::CA,start_index+3);

		calculate_strand_axis_long(c1,c2,c3,c4,startVec,endVec);

		return true;

	}

	PRODART::UTILS::vector3d_vector coords;

	coords.reserve(end_index - start_index);

	vector3d c1 = protein->get_bb_atom_coords(POSE::CA,start_index);
	vector3d c2 = protein->get_bb_atom_coords(POSE::CA,start_index+1);
	vector3d c3 = protein->get_bb_atom_coords(POSE::CA,start_index+2);
	vector3d c4 = protein->get_bb_atom_coords(POSE::CA,start_index+3);

	for (int i = start_index+4; i <= end_index; i++){

		vector3d  a,b;

		calculate_strand_axis(c1,c2,c3,c4,a,b);
		coords.push_back(a);
		coords.push_back(b);

		c1 = c2;
		c2 = c3;
		c3 = c4;
		c4 = protein->get_bb_atom_coords(POSE::CA,i);


	}

	const vector3d first_atom = protein->get_bb_atom_coords(POSE::CA,start_index);
	const vector3d last_atom = protein->get_bb_atom_coords(POSE::CA,end_index);
	PRODART::UTILS::line_fit3d(coords,
			startVec,
			endVec,
			first_atom,
			last_atom);

	return true;
	*/

}


















//! full extent helix calculation
bool get_helix_end_points_long(const PRODART::UTILS::vector3d_vector& ca_coords,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec){
	const unsigned int len = ca_coords.size();
	if (len <= 3){
		return false;
	}
	else if (len <= 4){

		//protein->get_bb_atom_coords(POSE::CA, 0);
		/*
		vector3d c1 = protein->get_bb_atom_coords(POSE::CA,start_index);
		vector3d c2 = protein->get_bb_atom_coords(POSE::CA,start_index+1);
		vector3d c3 = protein->get_bb_atom_coords(POSE::CA,start_index+2);
		vector3d c4 = protein->get_bb_atom_coords(POSE::CA,start_index+3);
		*/

		calculate_helix_axis_long(ca_coords[0],ca_coords[1],ca_coords[2],ca_coords[3],startVec,endVec);

		return true;

	}

	PRODART::UTILS::vector3d_vector coords;

	coords.reserve(len);

	vector3d c1 = ca_coords[0];//protein->get_bb_atom_coords(POSE::CA,start_index);
	vector3d c2 = ca_coords[1];//protein->get_bb_atom_coords(POSE::CA,start_index+1);
	vector3d c3 = ca_coords[2];//protein->get_bb_atom_coords(POSE::CA,start_index+2);
	vector3d c4 = ca_coords[3];//protein->get_bb_atom_coords(POSE::CA,start_index+3);

	for (unsigned int i = 4; i < len; i++){

		vector3d  a,b;

		calculate_helix_axis(c1,c2,c3,c4,a,b);
		coords.push_back(a);
		coords.push_back(b);

		c1 = c2;
		c2 = c3;
		c3 = c4;
		c4 = ca_coords[i];//protein->get_bb_atom_coords(POSE::CA,i);


	}

	const vector3d first_atom = ca_coords.front();//protein->get_bb_atom_coords(POSE::CA,start_index);
	const vector3d last_atom = ca_coords.back();//protein->get_bb_atom_coords(POSE::CA,end_index);
	PRODART::UTILS::line_fit3d(coords,
			startVec,
			endVec,
			first_atom,
			last_atom);

	return true;
}

//! full extent strand calculation
bool get_strand_end_points_long(const PRODART::UTILS::vector3d_vector& ca_coords,
		PRODART::UTILS::vector3d& startVec,
		PRODART::UTILS::vector3d& endVec){
	const unsigned int len = ca_coords.size();
	if (len <= 3){
		return false;
	}
	else if (len <= 4){

		//protein->get_bb_atom_coords(POSE::CA, 0);
		/*
		vector3d c1 = protein->get_bb_atom_coords(POSE::CA,start_index);
		vector3d c2 = protein->get_bb_atom_coords(POSE::CA,start_index+1);
		vector3d c3 = protein->get_bb_atom_coords(POSE::CA,start_index+2);
		vector3d c4 = protein->get_bb_atom_coords(POSE::CA,start_index+3);
		*/

		calculate_strand_axis_long(ca_coords[0],ca_coords[1],ca_coords[2],ca_coords[3],startVec,endVec);

		return true;

	}

	PRODART::UTILS::vector3d_vector coords;

	coords.reserve(len);

	vector3d c1 = ca_coords[0];//protein->get_bb_atom_coords(POSE::CA,start_index);
	vector3d c2 = ca_coords[1];//protein->get_bb_atom_coords(POSE::CA,start_index+1);
	vector3d c3 = ca_coords[2];//protein->get_bb_atom_coords(POSE::CA,start_index+2);
	vector3d c4 = ca_coords[3];//protein->get_bb_atom_coords(POSE::CA,start_index+3);

	for (unsigned int i = 4; i < len; i++){

		vector3d  a,b;

		calculate_strand_axis(c1,c2,c3,c4,a,b);
		coords.push_back(a);
		coords.push_back(b);

		c1 = c2;
		c2 = c3;
		c3 = c4;
		c4 = ca_coords[i];//protein->get_bb_atom_coords(POSE::CA,i);


	}

	const vector3d first_atom = ca_coords.front();//protein->get_bb_atom_coords(POSE::CA,start_index);
	const vector3d last_atom = ca_coords.back();//protein->get_bb_atom_coords(POSE::CA,end_index);
	PRODART::UTILS::line_fit3d(coords,
			startVec,
			endVec,
			first_atom,
			last_atom);

	return true;
}





}
}
