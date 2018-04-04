//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite.
//
// Copyright (c) 2010, James T. MacDonald
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// Neither the name of the James T. MacDonald nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
// or visit the vector3d project home: http://code.google.com/p/vector3d/
//
#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>

#include <limits>


namespace PRODART {

namespace UTILS {

//! use this PI definition in all calculations
const double PI = std::acos(-1); //3.14159265358979323846;

void dummy_funct(double val);

inline double degrees_to_radians(const double degrees){
	return PI * (degrees / 180.0);
}

inline double radians_to_degrees(const double radians){
	return 180.0 * (radians / PI);
}

//! returns true if "dih" (in radians) is between "to" and "from" in dihedral space
inline bool is_dihedral_in_range(const double from, const double to, const double dih){
	const double tol = degrees_to_radians(1.0);
	if (to >= from){
		if (dih >= from && dih <= to){
			return true;
		}
		else {
			return false;
		}
	}
	else {
		if (dih >= from && dih <= (PI+tol)){
			return true;
		}
		else if (dih >= -(PI + tol) && dih <= to){
			return true;
		}
		else {
			return false;
		}
	}
	return false;
}

inline double get_dihedral_distance(const double dih1, const double dih2){

	const double diff = dih1 - dih2;

	if (fabs(diff) < PI){
		return diff;
	}
	else if (diff > 0){
		return (diff - (2.0 * PI));
	}
	else {
		return (diff + (2.0 * PI));
	}

}

class quaternion_rotate_params;

//! NOTE: all angles are in RADIANS
class vector3d {

	friend std::ostream &operator<<( std::ostream&,  const vector3d & );

	friend vector3d operator-(const vector3d& a);

	//! calculate vector subtraction
	friend vector3d operator-( const vector3d&,  const vector3d& );
	//! calculate vector addition
	friend vector3d operator+( const vector3d&,  const vector3d& );
	//! calculate vector cross product
	friend vector3d operator*( const vector3d&, const vector3d& );
	friend vector3d operator/( const vector3d&, double d );
	friend vector3d operator*( const vector3d&, double d );
	friend vector3d operator*( double d, const vector3d& );

	friend bool operator==( const vector3d&,  const vector3d& );
	friend bool operator!=( const vector3d&,  const vector3d& );

	friend void quaternion_rotate(vector3d& v, const double angle, const vector3d& axis);
	friend void quaternion_rotate(vector3d& v, const quaternion_rotate_params& params);

public:
	vector3d();
	vector3d( double, double, double );
	vector3d(const vector3d&);

	vector3d& operator=(const vector3d&);
	vector3d& operator+=(const vector3d&);
	vector3d& operator-=(const vector3d&);
	vector3d& operator*=(const double val);
	double& operator[] (const int index);
	const double& operator[] (const int index) const;

	//! get vector magnitude squared
	double mod_sq() const;
	//! get vector magnitude
	double mod() const;
	//! calculate vector dot product
	double dot( const vector3d& ) const;
    // euclidean distance
	double dist(const vector3d& ) const;
    // euclidean distance squared
	double dist_sq(const vector3d& ) const;
	
	//! returns myVector squared
	vector3d sq() const;


	//! two ways to access vector components - either through the references x,y,z or through the operator[]:
	double &x, &y, &z;

private:
	double r[3];



};

typedef std::vector<vector3d> vector3d_vector;

class rot_matrix {

	friend void matrix_rotate(vector3d& v, const rot_matrix& rot);

public:

	//! identity matrix
	rot_matrix();
	rot_matrix(const rot_matrix&);
	rot_matrix& operator=(const rot_matrix&);

	/*! for use with qcprot routine where rot[9] is in the format:
	 * rot[9]   -- the rotation matrix in the order of xx, xy, xz, yx, yy, yz, zx, zy, zz
	 *
	*/
	rot_matrix(const double * const rot);

	double get(const int i, const int j) const {
		return m[i][j];
	}

	void set(const int i, const int j, double val) {
		m[i][j] = val;
	}

private:
	double m[3][3];

};

class quaternion_rotate_params {

public:
	double a ;//=  std::cos(angle / 2);
	double sin_constant ;//=  std::sin(angle / 2);

	double axis_0 ;//=  sin_constant*axis.r[0];
	double axis_1 ;//=  sin_constant*axis.r[1];
	double axis_2 ;//=  sin_constant*axis.r[2];

	double t2 ;//=    a*axis_0;
	double t3 ;//=    a*axis_1;
	double t4 ;//=    a*axis_2;
	double t5 ;//=   -axis_0*axis_0;
	double t6 ;//=    axis_0*axis_1;
	double t7 ;//=    axis_0*axis_2;
	double t8 ;//=   -axis_1*axis_1;
	double t9 ;//=    axis_1*axis_2;
	double t10 ;//=  -axis_2*axis_2;
	quaternion_rotate_params(const double angle, const vector3d& axis){
		a = std::cos(angle / 2);
		sin_constant = std::sin(angle / 2);

		axis_0 = sin_constant*axis[0];
		axis_1 = sin_constant*axis[1];
		axis_2 = sin_constant*axis[2];

		t2 =   a*axis_0;
		t3 =   a*axis_1;
		t4 =   a*axis_2;
		t5 =  -axis_0*axis_0;
		t6 =   axis_0*axis_1;
		t7 =   axis_0*axis_2;
		t8 =  -axis_1*axis_1;
		t9 =   axis_1*axis_2;
		t10 = -axis_2*axis_2;
	}
};



double angle(const vector3d& p1,
		const vector3d& p2,
		const vector3d& p3);


double dihedral(const vector3d& p1,
		const vector3d& p2,
		const vector3d& p3,
		const vector3d& p4);

vector3d dihedralEnd(const vector3d& p1,
		const vector3d& p2,
		const vector3d& p3,
		const double d34,
        const double a234,
        const double dihedral);

vector3d dihedralEnd_d24_a324(const vector3d& p1,
		const vector3d& p2,
		const vector3d& p3,
		const double d24,
        const double a324,
        const double dihedral);

double angle(const vector3d& p1, const vector3d& p2, const vector3d& p3);

double angle(const vector3d& p1, const vector3d& p2, const vector3d& q1, const vector3d& q2);

void quaternion_rotate(vector3d& v, const double angle, const vector3d& axis);
void quaternion_rotate(vector3d& v, const quaternion_rotate_params& params);

//! apply rotation matrix, rot, to vector3d, v - NOTE: needs verification
void matrix_rotate(vector3d& v, const rot_matrix& rot);

bool lineIntersect(const vector3d& p1,const vector3d& p2,const vector3d& p3,const vector3d& p4,
		vector3d& pa,vector3d& pb,
		double &mua, double &mub,
		double & dih_ang);

bool segmentIntersect(const vector3d& p1,const vector3d& p2,const vector3d& p3,const vector3d& p4,
		vector3d& pa,vector3d& pb,
		double &mua, double &mub,
		double & dih_ang,
		double & p1_pa_pb_ang,
		double & pa_pb_p3_ang,
		double & dist_sq);


/*********************************************************************************************
 *********************************************************************************************
 * Start of function definitions
 *********************************************************************************************
 *********************************************************************************************/





inline vector3d::vector3d( void ) : x( this->r[0]), y( this->r[1] ), z( this->r[2] ), r{0,0,0}  {

//	x=0;
//	y=0;
//	z=0;

}

inline vector3d::vector3d( double a, double b, double c ) : x( this->r[0]), y( this->r[1] ), z( this->r[2]) , r{a,b,c} {

//	x=a;
//	y=b;
//	z=c;

}

inline vector3d::vector3d(const vector3d& vec2) : x( this->r[0]), y( this->r[1] ), z( this->r[2] ){

	if( this !=  &vec2 ){
		this->r[0] = vec2.r[0];
		this->r[1] = vec2.r[1];
		this->r[2] = vec2.r[2];
	}

}

inline vector3d& vector3d::operator=(const vector3d& vec2) {

	if( this !=  &vec2 ){
		this->r[0] = vec2.r[0];
		this->r[1] = vec2.r[1];
		this->r[2] = vec2.r[2];
	}

	return *this;
}

inline vector3d& vector3d::operator+=(const vector3d& vec2){
	this->r[0] += vec2.r[0];
	this->r[1] += vec2.r[1];
	this->r[2] += vec2.r[2];

	return *this;
}

inline vector3d& vector3d::operator-=(const vector3d& vec2){
	this->r[0] -= vec2.r[0];
	this->r[1] -= vec2.r[1];
	this->r[2] -= vec2.r[2];

	return *this;
}

inline vector3d& vector3d::operator*=(const double val){
	this->r[0] *= val;
	this->r[1] *= val;
	this->r[2] *= val;

	return *this;
}

inline const double& vector3d::operator[] (const int index) const{
	return r[index];
}

inline double& vector3d::operator[] (const int index){
	return r[index];
}

inline double vector3d::mod_sq( void ) const {
	//! returns modulus (magnitude) squared
	return (x*x) + (y*y) + (z*z) ;
}

inline double vector3d::mod( void ) const {
	//! returns modulus (magnitude)
	return sqrt( (x*x) + (y*y) + (z*z) );
}

inline double vector3d::dot( const vector3d& vec2 ) const
{

	return (x * vec2.x) +  (y * vec2.y) +  (z * vec2.z);

}

inline double vector3d::dist( const vector3d& vec2) const{
	//double dx=x-vec2.x, dy=y-vec2.y, dz=z-vec2.z;
	return (*this - vec2).mod();//(sqrt(dx*dx + dy*dy + dz*dz));
}

inline double vector3d::dist_sq( const vector3d& vec2) const{
	return (*this - vec2).mod_sq();
}

inline vector3d vector3d::sq() const{
	return vector3d(x*x, y*y, z*z);
}

inline double angle(const vector3d& p1, const vector3d& p2, const vector3d& p3){
	const vector3d a = p2-p1;
    const vector3d b = p3-p2;

    //a = a/a.mod();
    //b = b/b.mod();

    return PI - acos(a.dot(b)/(a.mod()*b.mod()));
}


inline vector3d operator-(const vector3d& vec){
	vector3d returnVec;

	returnVec.x = - (vec.x);
	returnVec.y = - (vec.y);
	returnVec.z = - (vec.z);

	return returnVec;
}

inline vector3d operator-( const vector3d &vec1, const vector3d &vec2 )
{
	vector3d returnVec;

	returnVec.x = (vec1.x) - (vec2.x);
	returnVec.y = (vec1.y) - (vec2.y);
	returnVec.z = (vec1.z) - (vec2.z);

	return returnVec;
}


inline vector3d operator+( const vector3d &vec1, const vector3d &vec2 )
{
	vector3d returnVec;

	returnVec.x = (vec1.x) + (vec2.x);
	returnVec.y = (vec1.y) + (vec2.y);
	returnVec.z = (vec1.z) + (vec2.z);

	return returnVec;
}

inline vector3d operator*( const vector3d &vec1, const vector3d &vec2 )
{
	vector3d returnVec;

	returnVec.x = ((vec1.y)*(vec2.z)) - ((vec2.y)*(vec1.z));
    returnVec.y = ((vec1.z)*(vec2.x)) - ((vec2.z)*(vec1.x));
    returnVec.z = ((vec1.x)*(vec2.y)) - ((vec2.x)*(vec1.y));

	return returnVec;
}

inline vector3d operator/( const vector3d &vec1, double d )
{
	vector3d returnVec;

	returnVec.x = (vec1.x)/d;
	returnVec.y = (vec1.y)/d;
	returnVec.z = (vec1.z)/d;

	return returnVec;
}

inline vector3d operator*( const vector3d &vec1, double d )
{
	vector3d returnVec;

	returnVec.x = (vec1.x)*d;
	returnVec.y = (vec1.y)*d;
	returnVec.z = (vec1.z)*d;

	return returnVec;
}

inline vector3d operator*( double d, const vector3d &vec1 )
{
	vector3d returnVec;

	returnVec.x = (vec1.x)*d;
	returnVec.y = (vec1.y)*d;
	returnVec.z = (vec1.z)*d;

	return returnVec;
}

inline double dihedral(const vector3d& p1, const vector3d& p2, const vector3d& p3, const vector3d& p4){

	const vector3d a = p2 - p1;
    const vector3d b = p3 - p2;
    const vector3d c = p4 - p3;

    const vector3d nab = a*b;      // normal to plane ab
    const vector3d nbc = b*c;      // normal to plane bc

    const double mprod = (nab.mod()) * (nbc.mod());

	const double top = nab.dot(nbc);					//bug fix section
	double angle;
	fabs(top) < mprod ? angle= acos((top)/mprod) : ( top > 0.0 ? angle= 0.0 : angle = PI);	//test that top/mprod is not more than 1 as acos would return NaN

    return (a.dot(nbc)) >= 0 ? angle : -angle;


}



inline vector3d dihedralEnd(const vector3d& p1, const vector3d& p2, const vector3d& p3,
		const double d34,
		const double a234,
		const double dihedral){

    const vector3d a = p2-p1;
    const vector3d b = p3-p2;
    const vector3d nab = a*b;     // normal to plane a-b
    const vector3d n3 = nab*b;    // normal to plane nab-b

    ////  c = p4-p3 ////

    // component of c in the direction of b
    const vector3d c1 = (b/b.mod())*(-d34*cos(a234));

    // component of c in the direction of nab
    const vector3d c2 = (nab/nab.mod())*(d34*sin(a234)*sin(dihedral));

    // component of c in the direction normal to plane nab-b
    const vector3d c3 = (n3/n3.mod())*(d34*sin(a234)*cos(dihedral));

    return p3 + c1 + c2 + c3;
}

//! NEED to verify this routine
inline vector3d dihedralEnd_d24_a324(const vector3d& p1,
		const vector3d& p2,
		const vector3d& p3,
		const double d24,
        const double a324,
        const double dihedral){

    const vector3d a = p2-p1;
    const vector3d b = p3-p2;
    const vector3d nab = a*b;     // normal to plane a-b
    const vector3d n3 = nab*b;    // normal to plane nab-b

    ////  c = p4-p2 ////

    // component of c in the direction of b
    const vector3d c1 = (b/b.mod())*(d24*cos(a324));

    // component of c in the direction of nab
    const vector3d c2 = (nab/nab.mod())*(d24*sin(a324)*sin(dihedral));

    // component of c in the direction normal to plane nab-b
    const vector3d c3 = (n3/n3.mod())*(d24*sin(a324)*cos(dihedral));

    return p2 + c1 + c2 + c3;
}

inline double angle(const vector3d& p1, const vector3d& p2, const vector3d& q1, const vector3d& q2){
	vector3d a = p2-p1;
    vector3d b = q2-q1;

	return acos(a.dot(b)/(a.mod()*b.mod()));
}

inline void quaternion_rotate(vector3d& v, const quaternion_rotate_params& params){

	const double temp_x = 2*( (params.t8 + params.t10)*v.r[0] + (params.t6 -  params.t4)*v.r[1] + (params.t3 + params.t7)*v.r[2] ) + v.r[0];
	const double temp_y = 2*( (params.t4 +  params.t6)*v.r[0] + (params.t5 + params.t10)*v.r[1] + (params.t9 - params.t2)*v.r[2] ) + v.r[1];
	const double temp_z = 2*( (params.t7 -  params.t3)*v.r[0] + (params.t2 +  params.t9)*v.r[1] + (params.t5 + params.t8)*v.r[2] ) + v.r[2];

	v.r[0] = temp_x;
	v.r[1] = temp_y;
	v.r[2] = temp_z;

}


inline void quaternion_rotate(vector3d& v, const double angle, const vector3d& axis){
	//! needs to be checked and verified!!!!!!
	const quaternion_rotate_params params(angle, axis);

	/*
	const double a = std::cos(angle / 2);
	const double sin_constant = std::sin(angle / 2);

	const double axis_0 = sin_constant*axis.r[0];
	const double axis_1 = sin_constant*axis.r[1];
	const double axis_2 = sin_constant*axis.r[2];

	const double t2 =   a*axis_0;
	const double t3 =   a*axis_1;
	const double t4 =   a*axis_2;
	const double t5 =  -axis_0*axis_0;
	const double t6 =   axis_0*axis_1;
	const double t7 =   axis_0*axis_2;
	const double t8 =  -axis_1*axis_1;
	const double t9 =   axis_1*axis_2;
	const double t10 = -axis_2*axis_2;


	const double temp_x = 2*( (params.t8 + params.t10)*v.r[0] + (params.t6 -  params.t4)*v.r[1] + (params.t3 + params.t7)*v.r[2] ) + v.r[0];
	const double temp_y = 2*( (params.t4 +  params.t6)*v.r[0] + (params.t5 + params.t10)*v.r[1] + (params.t9 - params.t2)*v.r[2] ) + v.r[1];
	const double temp_z = 2*( (params.t7 -  params.t3)*v.r[0] + (params.t2 +  params.t9)*v.r[1] + (params.t5 + params.t8)*v.r[2] ) + v.r[2];

	v.r[0] = temp_x;
	v.r[1] = temp_y;
	v.r[2] = temp_z;
	*/
	quaternion_rotate(v, params);

}

inline void matrix_rotate(vector3d& v, const rot_matrix& rot){
	const vector3d orig_v = v;

	v.x = (rot.m[0][0] * orig_v.x) + (rot.m[0][1] * orig_v.y) + (rot.m[0][2] * orig_v.z);
	v.y = (rot.m[1][0] * orig_v.x) + (rot.m[1][1] * orig_v.y) + (rot.m[1][2] * orig_v.z);
	v.z = (rot.m[2][0] * orig_v.x) + (rot.m[2][1] * orig_v.y) + (rot.m[2][2] * orig_v.z);

}

inline std::ostream &operator<<( std::ostream &output, const vector3d &outVec )
{
	output << outVec.x << "\t"
		   << outVec.y << "\t"
		   << outVec.z << "\t";
	return output;
}

inline bool operator==( const vector3d &vec1, const vector3d &vec2 )
{


	return (vec1.x == vec2.x) && (vec1.y == vec2.y) && (vec1.z == vec2.z);
}

inline bool operator!=( const vector3d& vec1,  const vector3d& vec2 ){
	return !(vec1 == vec2);
}



}
}


#endif /*MYVECTOR_H_*/
