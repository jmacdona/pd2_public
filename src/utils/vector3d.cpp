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
#include "vector3d.h"


using std::cout;
using std::cin;
using std::endl;

using std::pow;
using std::sqrt;
using std::acos;
using std::asin;

using std::ostream;
using std::istream;


namespace PRODART {

namespace UTILS {

/*
bool is_dihedral_in_range(const double from, const double to, const double ang){
	const double tol = degrees_to_radians(1.0);
	if (to >= from){
		if (ang >= from && ang <= to){
			return true;
		}
		else {
			return false;
		}
	}
	else {
		if (ang >= from && ang <= (PI+tol)){
			return true;
		}
		else if (ang >= -(PI + tol) && ang <= to){
			return true;
		}
		else {
			return false;
		}
	}
	return false;
}

	if (!is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), degrees_to_radians(175))){
		cerr << "is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), 175) failed\n";
	}
	if (!is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), degrees_to_radians(-175))){
		cerr << "is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), -175) failed\n";
	}
	if (!is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), degrees_to_radians(180))){
		cerr << "is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), 180) failed\n";
	}
	if (!is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), degrees_to_radians(-180))){
		cerr << "is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), -180) failed\n";
	}
	if (!is_dihedral_in_range(degrees_to_radians(-90), degrees_to_radians(90), degrees_to_radians(-0))){
		cerr << "is_dihedral_in_range(degrees_to_radians(-90), degrees_to_radians(90), -0) failed\n";
	}
	if (is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), degrees_to_radians(0))){
		cerr << "is_dihedral_in_range(degrees_to_radians(170), degrees_to_radians(-170), -0) failed\n";
	}
	if (!is_dihedral_in_range(degrees_to_radians(-170), degrees_to_radians(170), degrees_to_radians(0))){
		cerr << "is_dihedral_in_range(degrees_to_radians(-170), degrees_to_radians(170), 0) failed\n";
	}
*/

void dummy_funct(double val){
	//do nothing
}

rot_matrix::rot_matrix(){
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			this->m[i][j] = 0;
		}
	}

	m[0][0] = 1;
	m[1][1] = 1;
	m[2][2] = 1;
}

rot_matrix::rot_matrix(const rot_matrix& rot2){

	if( this !=  &rot2 ){
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				this->m[i][j] = rot2.m[i][j];
			}
		}
	}

}

rot_matrix& rot_matrix::operator=(const rot_matrix& rot2){
	if( this !=  &rot2 ){
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				this->m[i][j] = rot2.m[i][j];
			}
		}
	}

	return *this;
}

rot_matrix::rot_matrix(const double * const rot){
	/*! for use with qcprot routine where rot[9] is in the format:
	 * rot[9]   -- the rotation matrix in the order of xx, xy, xz, yx, yy, yz, zx, zy, zz
	 *
	*/

	this->m[0][0] = rot[0];
	this->m[0][1] = rot[1];
	this->m[0][2] = rot[2];
	this->m[1][0] = rot[3];
	this->m[1][1] = rot[4];
	this->m[1][2] = rot[5];
	this->m[2][0] = rot[6];
	this->m[2][1] = rot[7];
	this->m[2][2] = rot[8];


}

bool lineIntersect(const vector3d& p1,const vector3d& p2,const vector3d& p3,const vector3d& p4,
		vector3d& pa,vector3d& pb,
		double &mua, double &mub,
		double & dih_ang){

	const double epsilon = std::numeric_limits<double>::min();

	vector3d p13,p43,p21;
	double d1343,d4321,d1321,d4343,d2121;
	double numer,denom;

	p13 = p1 - p3;

	p43 = p4 - p3;

	if (fabs(p43.x)  < epsilon && fabs(p43.y)  < epsilon && fabs(p43.z)  < epsilon)
		return false;

	p21 = p2 - p1;

	if (fabs(p21.x)  < epsilon && fabs(p21.y)  < epsilon && fabs(p21.z)  < epsilon)
		return false;

	d1343 = p13.dot(p43);//
	d4321 = p43.dot(p21);//
	d1321 = p13.dot(p21);//
	d4343 = p43.dot(p43);//
	d2121 = p21.dot(p21);//

	denom = d2121 * d4343 - d4321 * d4321;
	if (fabs(denom) < epsilon)
		return false;
	numer = d1343 * d4321 - d1321 * d4343;

	mua = numer / denom;
	mub = (d1343 + d4321 * (mua)) / d4343;



	pa = p1 + (mua * p21);

	pb = p3 + (mub * p43);

	dih_ang = dihedral(p1 + ((mua - 1.0) * p21), pa, pb, p3 + ((mub - 1.0) * p43));

	return true;

}


bool segmentIntersect(const vector3d& p1,const vector3d& p2,const vector3d& p3,const vector3d& p4,
		vector3d& pa,vector3d& pb,
		double &mua, double &mub,
		double & dih_ang,
		double & p1_pa_pb_ang,
		double & pa_pb_p3_ang,
		double & dist_sq){

	const double epsilon = std::numeric_limits<double>::min();

	vector3d   u = p2 - p1;//
	vector3d   v = p4 - p3;//
	vector3d   w = p1 - p3;//
	double    a = u.dot(u);//
	double    b = u.dot(v);//
	double    c = v.dot(v);//
	double    d = u.dot(w);//
	double    e = v.dot(w);//
	double    D = a*c - b*b;       // 
	double    sN, sD = D;      // 
	double    tN, tD = D;      // 


	if (D < epsilon) { // 
		sN = 0.0;        // 
		sD = 1.0;        // 
		tN = e;
		tD = c;
	}
	else {                //
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < 0.0) {       // 
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else if (sN > sD) {  // 
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN < 0.0) {           //
		tN = 0.0;
		// 
		if (-d < 0.0)
			sN = 0.0;
		else if (-d > a)
			sN = sD;
		else {
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD) {      //
		tN = tD;
		// 
		if ((-d + b) < 0.0)
			sN = 0;
		else if ((-d + b) > a)
			sN = sD;
		else {
			sN = (-d + b);
			sD = a;
		}
	}

	mua = (fabs(sN) < epsilon ? 0.0 : sN / sD);
	mub = (fabs(tN) < epsilon ? 0.0 : tN / tD);

	pa = p1 + (mua * u);

	pb = p3 + (mub * v);

	dist_sq = (pa-pb).mod_sq();

	dih_ang = dihedral(p1 + ((mua - 1.0) * u), pa, pb, p3 + ((mub - 1.0) * v));
	if (dist_sq > epsilon){
		p1_pa_pb_ang = angle(p1 + ((mua - 1.0) * u), pa, pb);
		pa_pb_p3_ang = angle(pa, pb, p3 + ((mub - 1.0) * v));
	}
	else {
		p1_pa_pb_ang = degrees_to_radians(90.0);
		pa_pb_p3_ang = degrees_to_radians(90.0);
	}



	return true;   //

}

}

}

/*
 *
 * most function definitions moved to header file so that they can be inlined
 *
 */
