//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
#include "line_fit.h"

using std::cout;
using std::cerr;
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

	bool pairCompareFirstGreater(const std::pair<double, int>& firstElem, const std::pair<double, int>& secondElem) {
		return firstElem.first > secondElem.first;
	}

bool is_vector3d_vector_valid(const vector3d_vector& vec){

	bool isOK = true;
	//vector3d_vector::const_iterator it;
	for (unsigned long it = 0; it < vec.size(); it++){
		const vector3d v3d = vec[it];
		for (int i = 0; i < 3; i++){
			if (boost::math::isinf(v3d[i])){
				cerr << "is_vector3d_vector_valid: ERROR: isinf at " << it << "[" << i << "]" << endl;
				isOK = false;
			}
			else if (boost::math::isnan(v3d[i])){
				cerr << "is_vector3d_vector_valid: ERROR: isnan at " << it << "[" << i << "]" << endl;
				isOK = false;
			}
		}
	}
	return isOK;
}
	
	vector3d get_CoM(const vector3d_vector& points){
		const unsigned long n = points.size();
		
		vector3d CoM(0,0,0);
		for (int i = 0; i < n; i++){
			CoM += points[i];
		}
		
		CoM = CoM / static_cast<double>(n);

		return CoM;
	}
	
	Eigen::Matrix3d get_A_matrix(const vector3d_vector& points, const vector3d CoM){
		const unsigned long n = points.size();
		Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
		for (int i = 0; i < n; i++){
			const vector3d c = points[i] - CoM;
			
			for (int j = 0 ; j < 3; j++){
				for (int k = j ; k < 3; k++){
					const double new_val = A(j,k) + (c[j] * c[k]);
					A(j,k) = new_val;
					if (j != k){
						A(k,j) = new_val;
					}
				}
			}
		}
		return A;
	}

	

	
	void eigen3_line_fit3d(const vector3d_vector& points,
						   vector3d& a,
						   vector3d& b){
		
		const unsigned long n = points.size();
		eigen3_line_fit3d(points, a, b, points[0], points[n-1]);
		
	}
	
	void eigen3_line_fit3d(const vector3d_vector& points,
						   vector3d& a,
						   vector3d& b,
						   const vector3d alternative_first,
						   const vector3d alternative_last){
		//const unsigned long n = points.size();
		
		vector3d CoM = get_CoM(points);
		Eigen::Matrix3d A = get_A_matrix(points, CoM);
		//cout << "Here is the matrix A:\n" << A << endl;
		
		
		
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(A);
		if (eigensolver.info() != Eigen::Success) {
			cerr << "ERROR: Eigen::SelfAdjointEigenSolver failed - eigensolver.info() != Eigen::Success";
		}
		/*
		cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << endl;
		cout << "Here's a matrix whose columns are eigenvectors of A \n"
		<< "corresponding to these eigenvalues:\n"
		<< eigensolver.eigenvectors() << endl;
		*/
		
		std::vector<std::pair<double, int> > sort_vec;
		for (int i = 0 ; i < 3; i++){
			sort_vec.push_back(std::pair<double, int>(eigensolver.eigenvalues()(i),i));
		}
		std::sort(sort_vec.begin(), sort_vec.end(), pairCompareFirstGreater);
		
		/*
		cout << "First Eigenval:\t" << sort_vec[0].first << endl;
		cout << "First Eigenvec:\t" << eigensolver.eigenvectors().col(sort_vec[0].second) << endl;
		*/
		
		vector3d v(eigensolver.eigenvectors().col(sort_vec[0].second)(0,0),
				   eigensolver.eigenvectors().col(sort_vec[0].second)(1,0),
				   eigensolver.eigenvectors().col(sort_vec[0].second)(2,0));
		
		v = v / v.mod();
		
		const vector3d first = alternative_first - CoM;
		const vector3d last = alternative_last - CoM;
		
		a = CoM + (first.dot(v) * v);
		b = CoM + (last.dot(v) * v);
		
		
	}
	
	void eigen3_line_fit3d(const vector3d_vector& points, vector3d& a, vector3d& b, double& dens, double& dvar, double* ev){
		
		vector3d CoM = get_CoM(points);
		Eigen::Matrix3d A = get_A_matrix(points, CoM);
		
		
		
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(A);
		if (eigensolver.info() != Eigen::Success) {
			cerr << "ERROR: Eigen::SelfAdjointEigenSolver failed - eigensolver.info() != Eigen::Success";
		}

		
		std::vector<std::pair<double, int> > sort_vec;
		for (int i = 0 ; i < 3; i++){
			sort_vec.push_back(std::pair<double, int>(eigensolver.eigenvalues()(i),i));
		}
		std::sort(sort_vec.begin(), sort_vec.end(), pairCompareFirstGreater);
		
		
		vector3d v(eigensolver.eigenvectors().col(sort_vec[0].second)(0,0),
				   eigensolver.eigenvectors().col(sort_vec[0].second)(1,0),
				   eigensolver.eigenvectors().col(sort_vec[0].second)(2,0));
		
		v = v / v.mod();
		
		const vector3d first = points[0] - CoM;
		const vector3d last = points[points.size()-1] - CoM;
		
		a = CoM + (first.dot(v) * v);
		b = CoM + (last.dot(v) * v);
		
		// find mean axial rise
		//
		double ds=0.0, e=0.0, lastd=0.0;
		vector3d c0=first.dot(v)*v;
		for(unsigned int k=1;k<points.size();k++)
		{
			vector3d pt=points[k]-CoM;
			vector3d p1=(pt.dot(v)*v);
			double d=p1.dist(c0);
			double f=d-lastd;
			//cout << "k: " << k << " c0: " << c0.x << " " << c0.y << " " << c0.z << " p1: " << p1.x << " " << p1.y << " " <<p1.z <<" D: " << f <<endl;
			ds+=f;
			e+=f*f;
			lastd=d;
		}
		dens=ds/static_cast<double>(points.size()-1);
		dvar=(e/static_cast<double>(points.size()-1)) - dens*dens;
		
		ev[0]=sort_vec[0].first;
		ev[1]=sort_vec[1].first;
		ev[2]=sort_vec[2].first;
		
	}


}
}




