/*
 * math_utils.h
 *
 *  Created on: May 28, 2012
 *      Author: jmacdona
 */

#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>

#include <limits>


namespace PRODART {
namespace UTILS {


//! increment vector when index 0 is least sig fig and index vec.size()-1 is most sig fig
template <class T>
inline bool incr_vector(std::vector<T>& vec, const std::vector<T>& upper_limits){
	const int num_ele = vec.size();
	//const int pose_len = rescount;

	for (int i = 0 ; i < num_ele; i++){
		if (vec[i] < upper_limits[i] -1){
			vec[i]++;
			break;
		}
		else if (i < num_ele-1){
			//reset if not last position
			vec[i] = 0;
		}
		else {
			// is last position and overflow
			vec[i]++;
			return false;
		}

	}

	bool overflow=true;
	for (int i = 0 ; i < num_ele; i++){
		overflow = overflow && (vec[i] >= upper_limits[i]);
	}
	if (overflow) return false;


	return true;
}



inline int int_pow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = int_pow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}




}
}

#endif /* MATH_UTILS_H_ */
