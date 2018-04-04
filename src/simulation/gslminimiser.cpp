//
//  gslminimiser.cpp
//  pd2
//
//  Created by jmacdon on 28/12/2014.
//  Copyright (c) 2014 James. All rights reserved.
//

#include <stdio.h>
#ifndef PD2_USE_CERES

#include "gslminimiser.h"





using namespace PRODART;
using namespace PRODART::UTILS;
using namespace PRODART::POSE;

using namespace std;


namespace PRODART {
	namespace POSE {
		namespace SIM {
			
			
			
			minimiser_interface_shared_ptr new_gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use){
				minimiser_interface_shared_ptr ptr(new gslminimiser(potentials_to_use));
				return ptr;
			}
			
			minimiser_interface_shared_ptr new_gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection){
				minimiser_interface_shared_ptr ptr(new gslminimiser(potentials_to_use, atom_selection));
				return ptr;
			}
			
			void gslminimiser::init(){
				lin_min_tol_ = ENV::get_option_value<double>("minimiser:lin_min_tol");//0.01;
				tol_ = ENV::get_option_value<double>("minimiser:tol");//1e-3;
				step_size = ENV::get_option_value<double>("minimiser:step_size");//0.01;
				
				if (tol_ < lin_min_tol_){
					cerr << "\nminimiser: WARNING!!!!: tol < lin_min_tol ---> " << tol_ << " < " <<  lin_min_tol_ << "\n" << endl;
				}
				
			}
			
			gslminimiser::gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector mask){
				init();
				this->potentials_ = potentials_to_use;
				b_min_selection_is_set = true;
				atom_selection_mask = mask;
				steps = 500;
				energy_map = potentials_->get_default_energies_map();
				
			}
			
			gslminimiser::gslminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use){
				init();
				this->potentials_ = potentials_to_use;
				b_min_selection_is_set = false;
				atom_selection_mask.resize(0);
				steps = 500;
				energy_map = potentials_->get_default_energies_map();
				
			}
			
			
			bool gslminimiser::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{
				this->pose_meta_ = meta_data;
				if (!initial_set_up()){
					return false;
				}
				
				return this->minimise();
			}
			
			bool gslminimiser::minimise() const{
				
				//ofstream output_pdb_file("traj.pdb", ios::out);
				
				const int pose_size = this->pose_meta_->get_pose()->get_all_atom_count();
				unsigned int x_size = pose_size * 3;
				
				
				gsl_vector* x = gsl_vector_alloc (x_size);
				//hold_vector = gsl_vector_alloc (x_size);
				
				assert(this->pose_meta_->get_pose()->coords_numerically_ok());
				
				this->pose_to_gsl_vector(this->pose_meta_->get_pose(), x);
				
				unsigned int iter = 0;
    int status = GSL_CONTINUE;
				
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
				
    gsl_multimin_function_fdf my_func;
				
    my_func.n = x_size;
    my_func.f = gslminimiser::f;
    my_func.df = gslminimiser::df;
    my_func.fdf = gslminimiser::fdf;
    my_func.params = (void*)this;
				
				// NOTE: 27/12/14: changing from gsl_multimin_fdfminimizer_vector_bfgs to gsl_multimin_fdfminimizer_vector_bfgs2
				// seems to substatially improve run times but this may only be
				// because it terminates sooner although it does seem to minimise in fewer steps...
    T = gsl_multimin_fdfminimizer_vector_bfgs;
    //T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, x_size);
				
    //gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);
				gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, lin_min_tol_);
				
    //energy_map.print_headers(cout);
				
				if (steps > 0){
					do
					{
						iter++;
						
						status = gsl_multimin_fdfminimizer_iterate (s);
						
						/*
						 cout << iter << "\t"
						 << s->f << "\t"
						 << endl;
						 */
						//energy_map.print_weighted_components(cout);
						
						if (status)
							break;
						
						status = gsl_multimin_test_gradient (s->gradient, tol_);
						
						
						
						//CA_mode_update_pdb(s->x);
						gsl_vector_to_pose(s->x, pose_meta_->get_pose());
						//pose_meta_->get_pose()->outputPdb(output_pdb_file);
						
						
						
					}
					while (status == GSL_CONTINUE && iter < steps);
				}
				
				cout << "min_steps_completed\t" << iter
				<< "\tmin_steps_requested\t" << steps << endl;
				
				/*
				 cout << iter << "\t"
				 << s->f << "\t"
				 << endl;
				 */
				
				gsl_multimin_fdfminimizer_free (s);
				
    gsl_vector_free (x);
				// gsl_vector_free (hold_vector);
				
    /*
	 energy_map.print_headers(cout);
	 energy_map.print_unweighted_components(cout);
	 energy_map.print_weights(cout);
	 energy_map.print_weighted_components(cout);
	 */
				
				return true;
			}
			
			bool  gslminimiser::initial_set_up() const{
				
				if (this->b_min_selection_is_set){
					const int sel_size = this->atom_selection_mask.size();
					const int pose_size = this->pose_meta_->get_pose()->get_all_atom_count();
					if (sel_size != pose_size){
						cerr << "ERROR: minimiser: atom_selection_mask difference size to pose." << endl;
						return false;
					}
					return true;
				}
				else {
					atom_selection_mask.resize(0);
					const int pose_size = this->pose_meta_->get_pose()->get_all_atom_count();
					atom_selection_mask.resize(pose_size, true);
					return true;
				}
				
			}
			
			
			double gslminimiser::f( const gsl_vector *x, void* obj ){
				gslminimiser& this_min = *((gslminimiser*)obj);
				this_min.gsl_vector_to_pose(x, this_min.pose_meta_->get_pose());
				this_min.pose_meta_->recalc_pair_lists_dists();
				const double energy =  this_min.potentials_->get_energy((this_min.pose_meta_), (this_min.energy_map));
				//cout << "f: " << energy  << endl;
				return energy;
			}
			void gslminimiser::df( const gsl_vector *x, void* obj, gsl_vector * g ){
				gslminimiser& this_min = *((gslminimiser*)obj);
				this_min.gsl_vector_to_pose(x, this_min.pose_meta_->get_pose());
				this_min.pose_meta_->recalc_pair_lists_dists();
				//const double energy =
				this_min.potentials_->get_energy_with_gradient((this_min.pose_meta_), (this_min.energy_map));
				this_min.vector3d_vector_to_gsl_vector(this_min.pose_meta_->get_gradient(), g);
				//cout << "df: " << energy  << endl;
			}
			void gslminimiser::fdf( const gsl_vector *x, void* obj, double *f, gsl_vector * g ){
				gslminimiser& this_min = *((gslminimiser*)obj);
				this_min.gsl_vector_to_pose(x, this_min.pose_meta_->get_pose());
				//this_min.pose_meta_->get_pose()->outputPdb(cout);
				this_min.pose_meta_->recalc_pair_lists_dists();
				*f = this_min.potentials_->get_energy_with_gradient((this_min.pose_meta_), (this_min.energy_map));
				this_min.vector3d_vector_to_gsl_vector(this_min.pose_meta_->get_gradient(), g);
				//cout << "fdf: " << *f  << endl;
			}
			
			void gslminimiser::pose_to_gsl_vector(PRODART::POSE::const_pose_shared_ptr pose_,  gsl_vector* coors) const{
				
				const int atom_count = pose_->get_all_atom_count();
				int gsl_vec_pos = 0;
				for (int i = 0 ; i < atom_count; i++){
					const vector3d coord = pose_->get_atom(i)->get_coords();
					for (int r = 0; r < 3; r++){
						/*
						 if(isnan(coord[r]) || isinf(coord[r])){
						 cerr << "minimiser: pose_to_gsl_vector: ERROR: is not a valid number: " << coord[r] << endl;
						 }
						 */
						gsl_vector_set(coors, gsl_vec_pos, coord[r]);
						gsl_vec_pos++;
					}
				}
				
			}
			void gslminimiser::gsl_vector_to_pose(const gsl_vector* coors, PRODART::POSE::pose_shared_ptr pose_) const{
				const int gsl_vec_size = coors->size;
				
				for (int i = 0 ; i < gsl_vec_size; i++){
					const int atom_num = i / 3;
					atom_shared_ptr atm = pose_->get_atom(atom_num);
					const int r = i % 3;
					const double val = gsl_vector_get(coors,i);
					/*
					 if(isnan(val) || isinf(val)){
					 cerr << "minimiser: gsl_vector_to_pose: ERROR: is not a valid number: " << val << endl;
					 }
					 */
					atm->set_r(r, val);
				}
				
				
			}
			void gslminimiser::vector3d_vector_to_gsl_vector(const PRODART::UTILS::vector3d_vector &grad,  gsl_vector* gslvec) const{
				const int atom_count = grad.size();
				int gsl_vec_pos = 0;
				for (int i = 0 ; i < atom_count; i++){
					if (this->atom_selection_mask[i] == true){
						const vector3d coord = grad[i];
						for (int r = 0; r < 3; r++){
							/*
							 if(isnan(coord[r]) || isinf(coord[r])){
							 cerr << "minimiser: vector3d_vector_to_gsl_vector: ERROR: is not a valid number: " << coord[r] << endl;
							 }
							 */
							gsl_vector_set(gslvec, gsl_vec_pos, coord[r]);
							gsl_vec_pos++;
						}
					}
					else {
						for (int r = 0; r < 3; r++){
							gsl_vector_set(gslvec, gsl_vec_pos, 0.0);
							gsl_vec_pos++;
						}
					}
				}
			}
			
			
			void gslminimiser::set_tol(double tol){
				tol_ = tol;
				
				if (tol_ < lin_min_tol_){
					lin_min_tol_ = tol_;
				}
				
			}
			
			void gslminimiser::set_lin_min_tol(double tol){
				lin_min_tol_ = tol;
			}
			
		}
	}
}

#endif 

