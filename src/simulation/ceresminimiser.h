//
//  ceresminimiser.h
//  pd2
//
//  Created by jmacdon on 28/12/2014.
//  Copyright (c) 2014 James. All rights reserved.
//

#ifndef __pd2__ceresminimiser__
#define __pd2__ceresminimiser__
#ifdef PD2_USE_CERES


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

#include "ceres/ceres.h"

#include "minimiser_interface.h"


namespace PRODART {
	namespace POSE {
		namespace SIM {
			
			typedef std::vector<bool> bool_vector;
			
			class ceresminimiser;
			
			typedef boost::shared_ptr<ceresminimiser> ceresminimiser_shared_ptr;
			
			minimiser_interface_shared_ptr new_ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use);
			minimiser_interface_shared_ptr new_ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection);
			
			//! gradient descent minimiser using ceres-solver library
			class ceresminimiser : public minimiser_interface {
				
				friend minimiser_interface_shared_ptr new_ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use);
				friend minimiser_interface_shared_ptr new_ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection);
				
			private:
				
				class minfunctor : public ceres::FirstOrderFunction {
				public:
					virtual ~minfunctor();
					minfunctor(PRODART::POSE::META::pose_meta_shared_ptr meta_data_,
							   PRODART::POSE::POTENTIALS::potentials_energies_map energy_map_,
							   PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use_,
							   bool_vector atom_selection_mask_);

					
					virtual bool Evaluate(const double* parameters,
										  double* cost,
										  double* gradient) const;
					
					virtual int NumParameters() const;
					

					bool param_to_pose(const double* parameters) const;
					bool pose_to_param( double* parameters) const;
					bool gradvec_to_gradient( const PRODART::UTILS::vector3d_vector grad_vec, double* gradient) const;

					
				private:
					typedef boost::tuple<unsigned int, unsigned int, POSE::atom_shared_ptr > atom_coord_tuple;
					typedef std::map<unsigned int, atom_coord_tuple> parameter_atom_coord_map;
					
					PRODART::POSE::META::pose_meta_shared_ptr meta_data;
					mutable PRODART::POSE::POTENTIALS::potentials_energies_map energy_map;
					PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use;
					bool_vector atom_selection_mask;
					parameter_atom_coord_map param_coord_map;
					
					bool setup_mapping();



				};
				
				
			protected:
				
				ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use);
				ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection);
				void init();
				
				bool initial_set_up() const;
				bool minimise() const;
				
				
				mutable PRODART::POSE::META::pose_meta_shared_ptr pose_meta_;
				PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_;
				bool b_min_selection_is_set;
				mutable bool_vector atom_selection_mask;
				
				unsigned int steps;
				double lin_min_tol_, tol_, step_size;
				
			public:
				
				bool make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const;
				
				void set_steps(unsigned int val){
					steps = val;
				}
				
				void set_tol(double tol);
				void set_lin_min_tol(double tol);
				
				
				
				PRODART::POSE::POTENTIALS::potentials_energies_map energy_map;
				
				
				
				
			};
			
			
			
			
			
		}
	}
}





#endif /* PD2_USE_CERES */
#endif /* defined(__pd2__ceresminimiser__) */
