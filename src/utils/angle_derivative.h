/*
 * angle_derivative.h
 *
 *  Created on: 1 Dec 2010
 *      Author: jmacdona
 */

#ifndef ANGLE_DERIVATIVE_H_
#define ANGLE_DERIVATIVE_H_
#include "vector3d.h"
#include <cassert>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>


namespace PRODART {
namespace UTILS {





inline void p1_theta_deriv(const vector3d & p1,
	const vector3d & p2,
	const vector3d & p3,
	vector3d & G){




	static const double min_angle( degrees_to_radians( double(0.1) ) );
	static const double max_angle( degrees_to_radians( double(179.9) ) );
	static const double max_x( std::cos( min_angle ));
	static const double min_x( std::cos( max_angle ));


	vector3d v1( p1 - p2 );
	vector3d v2( p3 - p2 );
	const double n1( v1.mod() );
	const double n2( v2.mod() );
	if ( n1 < double(1e-9) || n2 < double(1e-9) ) {
		return;
	}


	double x = (v1.dot(v2))/(n1*n2);



	vector3d f2(0.0, 0.0, 0.0);

	{
		const double f = double(1.0) / ( n1 * n2 );
		f2 += f * v2;
	}

	{
		const double f = double(-1.0) * x / ( n1 * n1 );
		f2 += f * v1;
	}

	x = std::min( std::max( min_x, x ), max_x );
	const double dtheta_dx = double(-1.0) / sqrt( double(1.0) - x*x );
	f2 *= dtheta_dx;



	assert(std::abs( ( f2.dot(v1) ) ) < double(1e-3));
	assert(std::abs( ( f2.dot( v1 * v2 ) ) ) < double(1e-3) );





	G += f2;
}







inline void angle_p1_deriv(const vector3d & p1,
		const vector3d & p2,
		const vector3d & p3,
		double & theta,
		vector3d & G){



	G = vector3d(0.0, 0.0, 0.0);
	theta = 0.0;

	vector3d u1( p1 - p2 );
	vector3d u2( p3 - p2 );
	const double n1_n2( u1.mod() * u2.mod() );
	if ( n1_n2 < double(1e-12) ) {
		std::cout << "WARNING: angle_p1_deriv: bond length: " << n1_n2 <<
				std::endl;
		std::cout << "WARNING: angle_p1_deriv:\tp1: " << p1
				<< "; p2: " << p2
				<< "; p3: " << p3 << std::endl;

		assert(!boost::math::isnan(n1_n2));
		assert(!boost::math::isinf(n1_n2));
		if (std::isnan(n1_n2) ||  std::isinf(n1_n2)){
			abort();
		}
		return;
	}

	p1_theta_deriv( p1, p2, p3,  G );

	double d( (u1.dot(u2)) / n1_n2 );
	const double tol(0.001);
	if ( d <= double(-1.0) + tol ) {

		d = double(-1.0) + tol;
	} else if ( d >= double(1.0) - tol ) {

		d = double(1.0) - tol;
	}
	theta = std::acos( d );

}




inline void angle_p2_deriv(const vector3d & p1,
		const vector3d & p2,
		const vector3d & p3,
		double & theta,
		vector3d & G){





	G = vector3d(0.0, 0.0, 0.0);
	theta = 0.0;

	vector3d v1( p1 - p2 );
	vector3d v2( p3 - p2 );
	const double v12( v1.mod() * v2.mod() );
	if ( v12 < double(1e-12) ) return;



	p1_theta_deriv( p2, p1, p3,  G );
	p1_theta_deriv( p2, p3, p1, G );


	assert( std::abs( ( G.dot( v1*v2) ) ) < double(1e-3) );

	double d( (v1.dot(v2)) / v12 );
	const double tol(0.001);
	if ( d <= double(-1.0) + tol ) {

		d = double(-1.0) +tol;
	} else if ( d >= double(1.0) - tol ) {

		d = double(1.0) -tol;
	}
	theta = std::acos( d );


	G *= double(-1.0);

}








}
}



#endif /* ANGLE_DERIVATIVE_H_ */
