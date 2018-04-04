//
//  bicubic_interp.h
//  bicubic_interp
//
//  Created by jmacdon on 10/01/2016.
//  Copyright (c) 2016 James MacDonald. All rights reserved.
//

#ifndef bicubic_interp_bicubic_interp_h
#define bicubic_interp_bicubic_interp_h

#include "boost/array.hpp"
#include "boost/multi_array.hpp"


class bicubic_interp {
	
public:
	typedef boost::multi_array<double, 2> mat_type;
	typedef boost::array<double, 4> array4_type;
	typedef boost::array<double, 16> array16_type;
	
	
	
public:
	mat_type c;
	double   x_lower,   x_upper,   y_lower,   y_upper;
	double x_size, y_size;
	
public:
	
	bicubic_interp(): c(boost::extents[4][4]),
	x_lower(0), x_upper(0), y_lower(0), y_upper(0),
	x_size(0), y_size(0){
		//c.assign(0);
	}
	
	// array4_type arrays are in the order 1:f(x_lower, y_lower), 2:f(x_upper, y_lower), 3:f(x_upper, y_upper), 4:f(x_lower, y_upper)
	bicubic_interp(array4_type f, array4_type df_dx, array4_type df_dy, array4_type d2f_dxdy,
				   double   x_lower_, double x_upper_, double y_lower_, double y_upper_) : c(boost::extents[4][4]),
	x_lower(x_lower_), x_upper(x_upper_), y_lower(y_lower_), y_upper(y_upper_),
	x_size(x_upper - x_lower), y_size(y_upper - y_lower){
		
		static int A_inv[16*16]=
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
			-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
			2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
			0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
			0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
			0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
			-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
			9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
			-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
			2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
			-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
			4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};
		

		const double area = x_size * y_size;
		array16_type cl, x;

		
		for (int i=0;i<4;i++) {
			x[i] = f[i];
			x[i+4] = df_dx[i]*x_size;
			x[i+8] = df_dy[i]*y_size;
			x[i+12] = d2f_dxdy[i]*area;
		}
		for (int i=0;i<16;i++) {
			double xx = 0;
			for (int k=0;k<16;k++) xx += A_inv[(16*i) + (k)] * x[k];
			cl[i] = xx;
		}
		
		int l = 0;
		for (int i=0; i<4; i++)
			for (int j=0; j<4; j++) c[i][j] = cl[l++];
		
	}
	
	double get_interp_val(const double x, const double y, double &df_dx, double &df_dy) const{

		double val = 0;
		df_dy = df_dx = 0.0;

		
		const double t = (x-x_lower)/x_size;
		const double u = (y-y_lower)/y_size;
		
		for (int i = 3;i >= 0;i--) {
			val  = t * val + ((c[i][3] * u + c[i][2]) * u + c[i][1]) * u +c[i][0];
			df_dy = t * df_dy + (3.0 * c[i][3] * u + 2.0 * c[i][2]) * u + c[i][1];
			df_dx = u * df_dx + (3.0 * c[3][i] * t + 2.0 * c[2][i]) * t + c[1][i];
		}
		
		df_dx /= x_size;
		df_dy /= y_size;
		
		return val;
	}
	
	
};

#endif
