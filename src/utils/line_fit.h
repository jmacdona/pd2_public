//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
#ifndef LINE_FIT_H_
#define LINE_FIT_H_

#include "vector3d.h"
#include <vector>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
#include <Eigen/Dense>



namespace PRODART {

namespace UTILS {

//! STL vector of vector3d objects
typedef std::vector< vector3d >  vector3d_vector;

bool is_vector3d_vector_valid(const vector3d_vector& vec);
	
vector3d get_CoM(const vector3d_vector& vec);
Eigen::Matrix3d get_A_matrix(const vector3d_vector& vec, const vector3d CoM);

	
//! do least squares fit to a bunch of points and return endpoints as vector3d a, b
	void eigen3_line_fit3d(const vector3d_vector& points,
					vector3d& return_first_endpoint,
					vector3d& return_last_endpoint);
	//! do least squares fit to a bunch of points and return endpoints as vector3d a, b
	inline void line_fit3d(const vector3d_vector& points,
					vector3d& return_first_endpoint,
						   vector3d& return_last_endpoint){
		eigen3_line_fit3d(points, return_first_endpoint, return_last_endpoint);
		return;
	}
	
	//! same as line_fit(const myVector_vector& , myVector& , myVector& ) but allows specification of different first, last points from endpoint projection
	void eigen3_line_fit3d(const vector3d_vector& points,
					vector3d& return_first_endpoint,
					vector3d& return_last_endpoint,
					const vector3d alternative_first,
					const vector3d alternative_last);
	//! same as line_fit(const myVector_vector& , myVector& , myVector& ) but allows specification of different first, last points from endpoint projection
	inline void line_fit3d(const vector3d_vector& points,
					vector3d& return_first_endpoint,
					vector3d& return_last_endpoint,
					const vector3d alternative_first,
					const vector3d alternative_last){
		eigen3_line_fit3d(points, return_first_endpoint, return_last_endpoint, alternative_first, alternative_last);
		return;
	}
	
	//! returns axial rise, variance around the axis and sorted first three eigenvalues
	void eigen3_line_fit3d(const vector3d_vector& points, vector3d& ep1, vector3d& ep2, double& dens, double& dvar, double* ev);
	//! returns axial rise, variance around the axis and sorted first three eigenvalues
	inline void line_fit3d(const vector3d_vector& points, vector3d& ep1, vector3d& ep2, double& dens, double& dvar, double* ev){
		eigen3_line_fit3d(points, ep1, ep2, dens, dvar, ev);
		return;
	}


}

}


#endif /*LINE_FIT_H_*/
