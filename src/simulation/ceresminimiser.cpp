//
//  ceresminimiser.cpp
//  pd2
//
//  Created by jmacdon on 28/12/2014.
//  Copyright (c) 2014 James. All rights reserved.
//

#include "ceresminimiser.h"
#ifdef PD2_USE_CERES





using namespace PRODART;
using namespace PRODART::UTILS;
using namespace PRODART::POSE;

using namespace std;


namespace PRODART {
	namespace POSE {
		namespace SIM {
			
			
			
			minimiser_interface_shared_ptr new_ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use){
				minimiser_interface_shared_ptr ptr(new ceresminimiser(potentials_to_use));
				return ptr;
			}
			
			minimiser_interface_shared_ptr new_ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector atom_selection){
				minimiser_interface_shared_ptr ptr(new ceresminimiser(potentials_to_use, atom_selection));
				return ptr;
			}
			
			void ceresminimiser::init(){
				lin_min_tol_ = ENV::get_option_value<double>("minimiser:lin_min_tol");//0.01;
				tol_ = ENV::get_option_value<double>("minimiser:tol");//1e-3;
				step_size = ENV::get_option_value<double>("minimiser:step_size");//0.01;
				steps = 500;

				if (tol_ < lin_min_tol_){
					cerr << "\nminimiser: WARNING!!!!: tol < lin_min_tol ---> " << tol_ << " < " <<  lin_min_tol_ << "\n" << endl;
				}
				
			}
			
			ceresminimiser::ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use, bool_vector mask){
				init();
				this->potentials_ = potentials_to_use;
				b_min_selection_is_set = true;
				atom_selection_mask = mask;
				energy_map = potentials_->get_default_energies_map();
				
			}
			
			ceresminimiser::ceresminimiser(PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use){
				init();
				this->potentials_ = potentials_to_use;
				b_min_selection_is_set = false;
				atom_selection_mask.resize(0);
				energy_map = potentials_->get_default_energies_map();
				
			}
			
			
			bool ceresminimiser::make_move(PRODART::POSE::META::pose_meta_shared_ptr meta_data) const{
				this->pose_meta_ = meta_data;
				if (!this->initial_set_up()){
					return false;
				}
				
				return this->minimise();
			}
			
			//actual minimisation happens here:
			bool ceresminimiser::minimise() const{
				
				cout << "USING NEW ceresminimiser!!!" << endl;
				
				//ofstream output_pdb_file("traj.pdb", ios::out);
				

				
				ceresminimiser::minfunctor *funct =new ceresminimiser::minfunctor(pose_meta_,
																			  energy_map,
																			  potentials_,
																			  atom_selection_mask);

				//const int pose_size = this->pose_meta_->get_pose()->get_all_atom_count();
				unsigned int param_size = funct->NumParameters();
				double *parameters = new double[param_size];//(double*) calloc (param_size,sizeof(double));
				//initialise parameters from pose
				funct->pose_to_param(parameters);
				
				ceres::GradientProblemSolver::Options options;
				options.minimizer_progress_to_stdout = false;//true;
				options.gradient_tolerance = this->tol_;
				options.min_line_search_step_size = step_size;
				options.max_num_iterations = steps;
				options.line_search_direction_type = ceres::LBFGS;
				options.function_tolerance = 1e-6; // 1e-6 is default
				
				ceres::GradientProblemSolver::Summary summary;
				
				//ceres::GradientProblem takes ownership of funct so assume it deletes afterwards
				ceres::GradientProblem problem(funct);
				ceres::Solve(options, problem, parameters, &summary);
				
				std::cout << summary.FullReport() << "\n";
				//std::cout << "Initial x: " << -1.2 << " y: " << 1.0 << "\n";
				//std::cout << "Final   x: " << parameters[0] << " y: " << parameters[1] << "\n";
				
				delete [] parameters;//free(parameters);
				
				return true;
			}
			
			bool  ceresminimiser::initial_set_up() const{
				
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
			
			



			
			void ceresminimiser::set_tol(double tol){
				tol_ = tol;
				
				if (tol_ < lin_min_tol_){
					lin_min_tol_ = tol_;
				}
				
			}
			
			void ceresminimiser::set_lin_min_tol(double tol){
				lin_min_tol_ = tol;
			}
			
			
			ceresminimiser::minfunctor::~minfunctor(){
			}
			
			ceresminimiser::minfunctor::minfunctor(PRODART::POSE::META::pose_meta_shared_ptr meta_data_,
												   PRODART::POSE::POTENTIALS::potentials_energies_map energy_map_,
												   PRODART::POSE::POTENTIALS::potentials_container_shared_ptr potentials_to_use_,
												   bool_vector atom_selection_mask_){
				this->meta_data = meta_data_;
				this->energy_map = energy_map_;
				this->potentials_to_use = potentials_to_use_;
				this->atom_selection_mask = atom_selection_mask_;
				this->setup_mapping();
			}
			
			bool ceresminimiser::minfunctor::Evaluate(const double* parameters,
													  double* cost,
													  double* gradient) const{
				//update pose coordinates from parameters
				this->param_to_pose(parameters);
				
				meta_data->recalc_pair_lists_dists();
				//const double energy =
				cost[0] = potentials_to_use->get_energy_with_gradient(meta_data, energy_map);
				gradvec_to_gradient(meta_data->get_gradient(), gradient);
				return true;
			}
			
			int ceresminimiser::minfunctor::NumParameters() const{
				//TODO change this:
				return (int)param_coord_map.size();
			}
			
			bool ceresminimiser::minfunctor::setup_mapping(){
				this->param_coord_map.clear();
				unsigned int ii;
				unsigned int param_index = 0;
				// ii is atom index
				for (ii = 0; ii < atom_selection_mask.size(); ii++){
					if (atom_selection_mask[ii] == true){
						atom_shared_ptr atm = meta_data->get_pose()->get_atom(ii);
						//x
						param_coord_map.insert(std::pair<unsigned int, boost::tuple<unsigned int, unsigned int, atom_shared_ptr> >(param_index, boost::tuple<unsigned int, unsigned int, atom_shared_ptr>(ii, 0, atm) ));
						//y
						param_coord_map.insert(std::pair<unsigned int, boost::tuple<unsigned int, unsigned int, atom_shared_ptr> >(param_index+1, boost::tuple<unsigned int, unsigned int, atom_shared_ptr>(ii, 1, atm) ));
						//z
						param_coord_map.insert(std::pair<unsigned int, boost::tuple<unsigned int, unsigned int, atom_shared_ptr> >(param_index+2, boost::tuple<unsigned int, unsigned int, atom_shared_ptr>(ii, 2, atm) ));
					
						param_index += 3;
					}
				}
				return true;
			}
			bool ceresminimiser::minfunctor::param_to_pose(const double* parameters) const{
				parameter_atom_coord_map::const_iterator it;
				for (it = param_coord_map.begin(); it != param_coord_map.end(); it++){
					unsigned int param_index = it->first;
					//unsigned int atom_index = it->second.get<0>();
					unsigned int coord_index = it->second.get<1>();
					atom_shared_ptr atm = it->second.get<2>();
					atm->set_r(coord_index, parameters[param_index]);
					//vector3d vec = atm->get_coords();
					//vec[coord_index] = parameters[param_index];
					//atm->set_coords(vec);
					
				}
				return true;
			}
			bool ceresminimiser::minfunctor::pose_to_param( double* parameters) const{
				parameter_atom_coord_map::const_iterator it;
				for (it = param_coord_map.begin(); it != param_coord_map.end(); it++){
					unsigned int param_index = it->first;
					//unsigned int atom_index = it->second.get<0>();
					unsigned int coord_index = it->second.get<1>();
					atom_shared_ptr atm = it->second.get<2>();
					parameters[param_index] = atm->get_r(coord_index);
					//vector3d vec = atm->get_coords();
					//vec[coord_index] = parameters[param_index];
					//atm->set_coords(vec);
					
				}
				return true;
			}
			
			bool ceresminimiser::minfunctor::gradvec_to_gradient( const PRODART::UTILS::vector3d_vector grad_vec, double* gradient) const{
				
				parameter_atom_coord_map::const_iterator it;
				for (it = param_coord_map.begin(); it != param_coord_map.end(); it++){
					unsigned int param_index = it->first;
					unsigned int atom_index = it->second.get<0>();
					unsigned int coord_index = it->second.get<1>();
					atom_shared_ptr atm = it->second.get<2>();
					gradient[param_index] = grad_vec[atom_index][coord_index];
					//vector3d vec = atm->get_coords();
					//vec[coord_index] = parameters[param_index];
					//atm->set_coords(vec);
					
				}
				
				return true;
			}
			
			
			
		}
	}
}








#endif /* PD2_USE_CERES */