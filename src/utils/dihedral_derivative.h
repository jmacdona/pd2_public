/*
 * dihedral_derivative.h
 *
 *  Created on: 6 Dec 2010
 *      Author: jmacdona
 */

#ifndef DIHEDRAL_DERIVATIVE_H_
#define DIHEDRAL_DERIVATIVE_H_
#include "vector3d.h"
#include <iostream>

namespace PRODART {
namespace UTILS {






inline void dih_helper(
	const vector3d  &v,
	const vector3d  &w,
	vector3d & G){

	const vector3d g = v * w;
	G += g;

}



inline void dihedral_deriv_second(
	const vector3d & p1,
	const vector3d & p2,
	const vector3d & p3,
	const vector3d & p4,
	double x,
	double & dih_ang,
	vector3d & G){


	const double min_angle( degrees_to_radians( double(0.1) ) );
	const double max_angle( degrees_to_radians( double(179.9) ) );
	const double max_x( std::cos( min_angle ));
	const double min_x( std::cos( max_angle ));




	dih_ang = dihedral( p1, p2, p3, p4 );

#ifndef NDEBUG

	const double dih_angU( std::acos( x ));

	const double dih_ang_diff =  ( std::abs( dih_ang ) - dih_angU );

	if ( !(dih_ang_diff < 1e-2) ){


		std::cerr << "dihedral_deriv_second: WARNING: std::abs( std::abs( dih_ang ) - dih_angU ) >= 1e-2: " <<  radians_to_degrees(dih_ang_diff) << "\t"
				<< "dih_ang: " << radians_to_degrees(dih_ang) << "\t"
				<< "dih_angU: " << radians_to_degrees(dih_angU) << "\t"
				<< "x: " << x << "\t"
				<< "p1: "
				<< p1 << "\tp2: "
				<< p2 << "\tp3: "
				<< p3 << "\tp4: "
				<< p4 << "\t"
				<< std::endl;

		//assert( std::abs( std::abs( dih_ang ) - dih_angU ) < 1e-2 );
	}
#endif



	x = std::min( std::max( min_x, x ), max_x );
	const double ddih_angU_dx = -1 / sqrt( 1- x*x );
	const double ddih_ang_ddih_angU( dih_ang < 0 ? -1 : 1 );

	G *= ddih_ang_ddih_angU * ddih_angU_dx;

}



inline void dihedral_p1_cosine_deriv_first(
	const vector3d & p1,
	const vector3d & p2,
	const vector3d & p3,
	const vector3d & p4,
	double & x,
	vector3d & G){


	G = vector3d(0.0, 0.0, 0.0);

	vector3d v1( p1-p2 );
	vector3d v2( p2-p3 );
	vector3d v3( p3-p4 );

	vector3d v12( ( v1 * v2 ));
	vector3d v23( ( v2 * v3 ));

	const double n12( v12.mod() );
	const double n23( v23.mod() );

	if ( n12 < double(1e-9) || n23 < double(1e-9) ) return;

	x = ( v12.dot( v23) ) / ( n12 * n23 );


	{
		const double f( double(1.0) / ( n12 * n23 ) );
		dih_helper( f * v2, v23 ,  G);
	}


	{
		const double f( double(-1.0) * x / ( n12 * n12 ) );
		dih_helper(  f * v2, v12,  G );
	}



	assert( std::abs( ( G.dot(v1) ) ) < double(1e-3) );
	assert( std::abs( ( G.dot(v2) ) ) < double(1e-3) );

}



inline void dihedral_p1_cosine_deriv(
	const vector3d & p1,
	const vector3d & p2,
	const vector3d & p3,
	const vector3d & p4,
	double & dih_ang,
	vector3d & G){


	double x( 0.0 );
	dihedral_p1_cosine_deriv_first( p1, p2, p3, p4, x, G );
	dihedral_deriv_second( p1, p2, p3, p4, x, dih_ang, G );

}



inline void dihedral_p2_cosine_deriv_first(
	const vector3d & p1,
	const vector3d & p2,
	const vector3d & p3,
	const vector3d & p4,
	double & x,
	vector3d & G){


	G = vector3d(0.0, 0.0, 0.0);


	vector3d v1( p1-p2 );
	vector3d v2( p2-p3 );
	vector3d v3( p3-p4 );

	vector3d v12( ( v1 * v2 ));
	vector3d v23( ( v2 * v3 ));

	const double n12( v12.mod() );
	const double n23( v23.mod() );

	if ( n12 < double(1e-9) || n23 < double(1e-9) ) return;

	x = ( v12.dot(v23)) / ( n12 * n23 );



	{

		{
			const double f( double(-1.0)/ ( n12 * n23 ) );
			dih_helper( f * v2, v23 , G);
		}

		{
			const double f( double(-1.0)/ ( n12 * n23 ) );
			dih_helper(  f * v1, v23 ,  G);
		}

		{
			const double f( double(1.0)/ ( n12 * n23 ) );
			dih_helper(  f * v3, v12 , G);
		}
	}

	{

		{
			const double f( x / ( n12 * n12 ) );
			dih_helper(  f * v2, v12,  G );
		}

		{
			const double f( x / ( n12 * n12 ) );
			dih_helper(  f * v1, v12,  G );
		}

		{
			const double f( double(-1.0) * x / ( n23 * n23 ) );
			dih_helper(  f * v3, v23, G );
		}
	}


	assert( std::abs( ( G.dot(v2) ) ) < double(1e-3) );

}



inline void dihedral_p2_cosine_deriv(
	const vector3d & p1,
	const vector3d & p2,
	const vector3d & p3,
	const vector3d & p4,
	double & dih_ang,
	vector3d & G){


	double x( 0.0 );
	dihedral_p2_cosine_deriv_first( p1, p2, p3, p4, x,  G );
	dihedral_deriv_second( p1, p2, p3, p4, x, dih_ang,  G );

}








}
}

#endif /* DIHEDRAL_DERIVATIVE_H_ */
