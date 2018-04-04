//
//  gslminimiser.h
//  pd2
//
//  Created by jmacdon on 28/12/2014.
//  Copyright (c) 2014 James. All rights reserved.
//

#ifndef pd2_gslminimiser_h
#define pd2_gslminimiser_h


#include "pose/pose.h"
#include "utils/line_fit.h"
#include "utils/vector3d.h"
#include "pose_utils/qcprot.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/pose_basic_kine.h"
#include "pose_meta/pose_meta_interface.h"
#include "pose_meta/ca_pose_meta.h"
#include "potentials/potentials_name.h"
#include "potentials/potentials_store.h"
#include "potentials/potential_interface.h"
#include "potentials/potentials_container.h"
#include "prodart_env/prodart_env.h"
#include "movers/mover_interface.h"
#include "movers/move_set.h"
#include <string>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>


#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "minimiser_interface.h"


namespace PRODART {
	namespace POSE {
		namespace SIM {
			
			typedef std::vector<bool> bool_vector;
			
			class gslminimiser;
			
			typedef boost::shared_ptr<gslminimiser> minimiser_shared_ptr;
			
			minimiser_interface_shared_ptr new_gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use);
			minimiser_interface_shared_ptr new_gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection);
			
			//! gradient descent minimiser using gsl library
			class gslminimiser : public minimiser_interface {
				
				friend minimiser_interface_shared_ptr new_gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use);
				friend minimiser_interface_shared_ptr new_gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection);
				
			private:
				
				
				
			protected:
				
				gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use);
				gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection);
				void init();
				
				bool initial_set_up() const;
				bool minimise() const;
				
				void pose_to_gsl_vector(PRODART::POSE::const_pose_shared_ptr pose_,  gsl_vector* coors) const;
				void gsl_vector_to_pose(const gsl_vector* coors, PRODART::POSE::pose_shared_ptr pose_) const;
				void vector3d_vector_to_gsl_vector(const PRODART::UTILS::vector3d_vector &grad,  gsl_vector* gslvec) const;
				
				mutable PRODART::POSE::META::pose_meta_shared_ptr pose_meta_;
				PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_;
				bool b_min_selection_is_set;
				mutable bool_vector atom_selection_mask;
				
				unsigned int steps;
				double lin_min_tol_,
				tol_,
				step_size;
				
			public:
				
				bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const;
				
				void set_steps(unsigned int val){
					steps = val;
				}
				
				void set_tol(double tol);
				void set_lin_min_tol(double tol);
				
				//static functions for use with GSL
				static double f( const gsl_vector *x, void* obj );
				static void df( const gsl_vector *x, void* obj, gsl_vector * g );
				static void fdf( const gsl_vector *x, void* obj, double *f, gsl_vector * g );
				
				
				
				PRODART::POSE::POTENTIALS::potentials_energies_map energy_map;
				
				
				
				
			};
			
			
			
			
			
		}
	}
}








#endif
