//
//  hist2d_interp.h
//  bicubic_interp
//
//  Created by jmacdon on 16/01/2016.
//  Copyright (c) 2016 James. All rights reserved.
//

#ifndef bicubic_interp_hist2d_interp_h
#define bicubic_interp_hist2d_interp_h

#include "hist_utils.h"
#include "bicubic_interp.h"
#include "multivariate_normal.h"

class hist2d_interp : public multi_hist<2>{

public:

	typedef boost::multi_array<bicubic_interp, 2> bicubic_interp_2darray_type;

	hist_binner_array_type bicubic_binners;
	bicubic_interp_2darray_type bicubic_arr;

	bool periodic_boundaries, always_recalc_interp;

	hist2d_interp() : multi_hist(), bicubic_binners(), bicubic_arr(), periodic_boundaries(false), always_recalc_interp(false)    {
	}

	hist2d_interp(const hist2d_interp& orig_) : multi_hist(orig_), bicubic_binners(orig_.bicubic_binners), periodic_boundaries(orig_.periodic_boundaries), always_recalc_interp(orig_.always_recalc_interp)    {
	    bicubic_arr.resize(boost::extents[num_bins[0]+1][num_bins[1]+1]);
	    bicubic_arr = orig_.bicubic_arr;
	}

	hist2d_interp& operator=(const hist2d_interp& orig_) {
	    multi_hist<2>::operator=(orig_);

	    bicubic_binners = orig_.bicubic_binners;
        periodic_boundaries = orig_.periodic_boundaries;
        always_recalc_interp = orig_.always_recalc_interp;

	    bicubic_arr.resize(boost::extents[num_bins[0]+1][num_bins[1]+1]);
	    bicubic_arr = orig_.bicubic_arr;

	    return *this;
	}


	hist2d_interp(uint_array_type num_bins_, double_array_type min_vals_, double_array_type max_vals_, double incr_, double pseudo_count,
			const bool periodic_ , const bool recalc_) :
		multi_hist(num_bins_, min_vals_, max_vals_, incr_, pseudo_count),
		bicubic_arr(boost::extents[num_bins_[0]+1][num_bins_[1]+1]),
		periodic_boundaries(periodic_), always_recalc_interp(recalc_){

		// set up binners for bicubic interp
		bicubic_binners[0] = hist_binner(binners[0].get_num_bins()+1,
				binners[0].get_min_val() - (binners[0].get_bin_size()/2),
				binners[0].get_max_val() + (binners[0].get_bin_size()/2.0));
		bicubic_binners[1] = hist_binner(binners[1].get_num_bins()+1,
				binners[1].get_min_val() - (binners[1].get_bin_size()/2),
				binners[1].get_max_val() + (binners[1].get_bin_size()/2.0));

	}


	uint_array_type get_bicubic_indicies(const double x, const double y) const{
		uint_array_type indicies;
		indicies[0] = bicubic_binners[0].get_bin(x);
		indicies[1] = bicubic_binners[1].get_bin(y);
		return indicies;

	}

	double get_interp_val(const double x, const double y, double &df_dx, double &df_dy){
		double_array_type pos;
		pos[0] = x;
		pos[1] = y;
		uint_array_type bicubic_indicies = get_bicubic_indicies(x, y);

		if (always_recalc_interp){
			bicubic_interp bi =  make_interp_obj(bicubic_indicies);
			const double val = bi.get_interp_val(x, y, df_dx, df_dy);
			return val;
		}
		else {
			const double val = bicubic_arr(bicubic_indicies).get_interp_val(x, y, df_dx, df_dy);
			return val;
		}

		return 0;
	}

	double get_interp_val_no_recalc(const double x, const double y, double &df_dx, double &df_dy) const{
			double_array_type pos;
			pos[0] = x;
			pos[1] = y;
			uint_array_type bicubic_indicies = get_bicubic_indicies(x, y);

			const double val = bicubic_arr(bicubic_indicies).get_interp_val(x, y, df_dx, df_dy);
			return val;

			return 0;
		}

	void precal_bicubic_arr(){
		for (unsigned int ii = 0; ii < bicubic_binners[0].get_num_bins(); ii++  ){
			for (unsigned int jj = 0; jj < bicubic_binners[0].get_num_bins(); jj++  ){
				uint_array_type indicies = {ii, jj};
				bicubic_arr(indicies) = make_interp_obj(indicies);
			}
		}
	}


	double get_f(const int i, const int j, const bool is_periodic){
		//std::cout << "DEB: calling get_f():\ti:" <<i << "\tj:" << j << std::endl;
        const int i_max = num_bins[0];
		const int j_max = num_bins[1];

		// periodic case
		if (is_periodic){
			return hist(make_indices( (i + i_max) %  i_max , (j + j_max) % j_max));
		}
		// non-periodic edge cases
		else {
			if (i >= 0 && j >= 0 && i < i_max && j < j_max){
			    return hist(make_indices(i,j));
			}
			else if (i == -1 && j >= 0 && j < j_max){
				const double diff = (get_f(1,j, is_periodic) - hist(make_indices(0,j)));
				const double f_val = hist(make_indices(0,j)) - diff;
				return f_val;
			}
			else if (i == i_max && j >= 0 && j < j_max){
				const double diff = get_f(i_max - 1, j, is_periodic) - get_f(i_max - 2, j, is_periodic);
				const double f_val = get_f(i_max - 1, j, is_periodic) + diff;
				return f_val;
			}
			else if (j == -1 && i >= 0 && i < i_max){
				const double diff = (get_f(i,1, is_periodic) - hist(make_indices(i,0)));
				const double f_val = hist(make_indices(i,0)) - diff;
				return f_val;
			}
			else if (j == j_max && i >= 0 && i < i_max){
				const double diff = (get_f(i, j_max - 1, is_periodic) - get_f(i, j_max - 2, is_periodic));
				const double f_val = hist(make_indices(i, j_max - 1)) + diff;
				return f_val;
			}
			else if (i >= -1 && j >= -1 && i <= i_max && j <= j_max) {
				// both i and j off the edge
				// gonna simply take the mean of neighbouring edge values - probably a better way of doing this
				const double f_i = get_f(std::min(std::max(i, 0), i_max - 1), j, is_periodic);
				const double f_j = get_f(i, std::min(std::max(j,0), j_max - 1), is_periodic);
				const double f_val = (f_i + f_j) / 2;
				//std::cout << "i:" <<i << "\tj:" << j << "\tf_i:" << f_i << "\tf_j:" << f_j << "\tf_val:" << f_val  << std::endl;
				return f_val;
			}
		}


		// catch all
		//std::cerr << "DEBUG: get_f(): reaching catch all condition: i:" << i << " j:" << j << std::endl;
		return hist(make_indices(std::min(std::max(i, 0), i_max - 1), std::min(std::max(j,0), j_max - 1)));


	}

	void get_gradients(const int_array_type indices, double &f, double &df_dx, double &df_dy, double &d2f_dxdy){
		//const int i_max = num_bins[0];
		//const int j_max = num_bins[1];
		const int i = indices[0];
		const int j = indices[1];
		int i_m1 = i-1, j_m1 = j-1, i_p1 = i+1, j_p1 = j+1;

		double val_i_p1_j=0, val_i_m1_j=0,
				val_i_j_p1=0, val_i_j_m1=0,
				val_i_p1_j_p1=0, val_i_p1_j_m1=0, val_i_m1_j_p1=0, val_i_m1_j_m1=0,
				x_width=0, y_width = 0;

		//if (i >= 0 && j >= 0 && i < i_max && j < j_max){
		f = get_f(i, j, periodic_boundaries);// hist(indices);
		//}



		val_i_p1_j = (get_f(i_p1, j, periodic_boundaries));
		val_i_m1_j = (get_f(i_m1, j, periodic_boundaries));
		val_i_j_p1 = (get_f(i, j_p1, periodic_boundaries));
		val_i_j_m1 = (get_f(i, j_m1, periodic_boundaries));

		val_i_p1_j_p1 = (get_f(i_p1, j_p1, periodic_boundaries));
		val_i_p1_j_m1 = (get_f(i_p1, j_m1, periodic_boundaries));
		val_i_m1_j_p1 = (get_f(i_m1, j_p1, periodic_boundaries));
		val_i_m1_j_m1 = (get_f(i_m1, j_m1, periodic_boundaries));

		x_width = binners[0].get_bin_size() * 2.0;
		y_width = binners[1].get_bin_size() * 2.0;

		df_dx = (val_i_p1_j - val_i_m1_j) / x_width;
		df_dy = (val_i_j_p1 - val_i_j_m1) / y_width;
		d2f_dxdy = (val_i_p1_j_p1 - val_i_p1_j_m1 - val_i_m1_j_p1 + val_i_m1_j_m1)
									/ ( x_width *   y_width );

	}

	bicubic_interp make_interp_obj(const uint_array_type bin_indices){

		const int i = bin_indices[0];
		const int j = bin_indices[1];

		bicubic_interp::array4_type f = {0,0,0,0},
				df_dx = {0,0,0,0},
				df_dy = {0,0,0,0},
				d2f_dxdy = {0,0,0,0};
		const double x_lower = bicubic_binners[0].get_bin_lowerbound(bin_indices[0]);
		const double x_upper = bicubic_binners[0].get_bin_upperbound(bin_indices[0]);
		const double y_lower = bicubic_binners[1].get_bin_lowerbound(bin_indices[1]);
		const double y_upper = bicubic_binners[1].get_bin_upperbound(bin_indices[1]);

		// check for standard case
		const int i_m1 = i -1, j_m1 = j - 1;
		const int_array_type indicies_0 = {(i_m1), (j_m1)},
				indicies_1 = {(i), (j_m1)},
				indicies_2 = {(i), (j)},
				indicies_3 = {(i_m1), (j)};

		double temp_f, temp_df_dx = 0, temp_df_dy = 0, temp_d2f_dxdy = 0;

		get_gradients(indicies_0, temp_f, temp_df_dx, temp_df_dy, temp_d2f_dxdy);
		f[0] = temp_f;
		df_dx[0] = temp_df_dx;
		df_dy[0] = temp_df_dy;
		d2f_dxdy[0] = temp_d2f_dxdy;

		get_gradients(indicies_1, temp_f, temp_df_dx, temp_df_dy, temp_d2f_dxdy);
		f[1] = temp_f;
		df_dx[1] = temp_df_dx;
		df_dy[1] = temp_df_dy;
		d2f_dxdy[1] = temp_d2f_dxdy;

		get_gradients(indicies_2, temp_f, temp_df_dx, temp_df_dy, temp_d2f_dxdy);
		f[2] = temp_f;
		df_dx[2] = temp_df_dx;
		df_dy[2] = temp_df_dy;
		d2f_dxdy[2] = temp_d2f_dxdy;

		get_gradients(indicies_3, temp_f, temp_df_dx, temp_df_dy, temp_d2f_dxdy);
		f[3] = temp_f;
		df_dx[3] = temp_df_dx;
		df_dy[3] = temp_df_dy;
		d2f_dxdy[3] = temp_d2f_dxdy;


		bicubic_interp bi = bicubic_interp(f, df_dx, df_dy, d2f_dxdy, x_lower, x_upper, y_lower, y_upper);


		/*
		std::cout << "DEB: making obj for i:" << i << " j:" << j << std::endl;
		std::cout << "DEB: x_bounds:\t" << bi.x_lower << "\t" << bi.x_upper << std::endl;
		std::cout << "DEB: y_bounds:\t" << bi.y_lower << "\t" << bi.y_upper << std::endl;
		std::cout << "DEB: f:\t";
		for (int ii = 0 ; ii < 4; ii++){
			std::cout << "\t" << f[ii];
		}
		std::cout << std::endl;
		std::cout << std::endl;
		*/


		return bi;
	}


    double calc_boundary_correction_term(const double val, const double stddev, const double L, const double U) const{
        // see J. Chem. Phys. 139:084102 (2013)
        //const double const_term = (std::sqrt(M_PI / 2));
        //const double var_term = stddev / (U - L);
        const double erf_term1 = std::erf((val - L)/(std::sqrt(2.0)*stddev));
        const double erf_term2 = std::erf((U - val)/(std::sqrt(2.0)*stddev));
        const double erf_term = erf_term1 + erf_term2;//std::erf((val - L)/(std::sqrt(2.0)*stddev)) + std::erf((U - val)/(std::sqrt(2.0)*stddev));
        //const double rtn_val_original = const_term * var_term * erf_term;
        const double rtn_val_scaled = 0.5 * erf_term;
        //std::cout << "corr: " << const_term << "\t" << var_term << "\t" << erf_term << "\t" << rtn_val_original << "\t" << rtn_val_scaled << std::endl;
        return rtn_val_scaled;//rtn_val_original;

    }

    void new_add_gaussian_count (const double weight, double_array_type p, double_array_type stddev, const double stddev_cutoff = 3.0, const bool fast = true){
        uint_array_type indicies = get_indicies(p);

        double corr_stddev_0 = stddev[0];//(std::min(std::min(stddev[0], p[0]), binners[0].get_max_val() - p[0]));
        double corr_stddev_1 = stddev[1];//(std::min(std::min(stddev[1], p[1]), binners[1].get_max_val() - p[1]));

        if (corr_stddev_0 <= (binners[0].get_bin_size()/stddev_cutoff) && corr_stddev_1 <= (binners[1].get_bin_size()/stddev_cutoff)){
            //std::cout << "adding count " << weight << " to: " << indicies[0] << " " << indicies[1] << " current_count: " << get_count(indicies) << std::endl;
            add_count(indicies, weight);
            //std::cout << "new count: " << get_count(indicies) << std::endl;
            return;
        }
        else
        {
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar = Eigen::MatrixXd::Zero(2, 2);
            Eigen::Matrix<double,Eigen::Dynamic, 1> mean = Eigen::MatrixXd::Zero(2, 1);

            covar(0,0) = corr_stddev_0 * corr_stddev_0;
            covar(1,1) = corr_stddev_1 * corr_stddev_1;
            mean(0) = p[0];
            mean(1) = p[1];

            MultivariateNormalPDF pdf;
            pdf.setMean(mean);
            pdf.setCovar(covar);

            const double cutoff_range_0 = stddev_cutoff * corr_stddev_0;
            const double cutoff_range_1 = stddev_cutoff * corr_stddev_1;

            const unsigned int cutoff_bin_range_0 = std::ceil(cutoff_range_0 / binners[0].get_bin_size());
            const unsigned int cutoff_bin_range_1 = std::ceil(cutoff_range_1 / binners[1].get_bin_size());

            const unsigned int ii_lower = (int(indicies[0]) - int(cutoff_bin_range_0)) > 0 ? indicies[0] - cutoff_bin_range_0 : 0;
            const unsigned int jj_lower = (int(indicies[1]) - int(cutoff_bin_range_1)) > 0 ? indicies[1] - cutoff_bin_range_1 : 0;

            const unsigned int ii_upper = (indicies[0] + cutoff_bin_range_0  ) < num_bins[0] ? indicies[0] + cutoff_bin_range_0  : num_bins[0]-1;
            const unsigned int jj_upper = (indicies[1] + cutoff_bin_range_1  ) < num_bins[1] ? indicies[1] + cutoff_bin_range_1  : num_bins[1]-1;


            Eigen::Matrix<double,Eigen::Dynamic, 1> pos = Eigen::MatrixXd::Zero(2, 1);
            Eigen::Matrix<double,Eigen::Dynamic, 1> pos_UU = Eigen::MatrixXd::Zero(2, 1);
            Eigen::Matrix<double,Eigen::Dynamic, 1> pos_UL = Eigen::MatrixXd::Zero(2, 1);
            Eigen::Matrix<double,Eigen::Dynamic, 1> pos_LU = Eigen::MatrixXd::Zero(2, 1);
            Eigen::Matrix<double,Eigen::Dynamic, 1> pos_LL = Eigen::MatrixXd::Zero(2, 1);



            /*
            double total_add = 0;
            for (unsigned int ii = ii_lower; ii <= ii_upper; ii++){
                for (unsigned int jj = jj_lower; jj <= jj_upper; jj++){
                    pos(0) = binners[0].get_bin_midpoint(ii);
                    pos(1) = binners[1].get_bin_midpoint(jj);
                    double val = this->bin_volume * pdf.getPdf(pos);
                    total_add+= weight*(val);
                    //add_count(make_indices(ii,jj), weight*(val));
                }
            }
            */



            const double norm_term = 1;//weight / total_add;


            for (unsigned int ii = ii_lower; ii <= ii_upper; ii++)
            {
                for (unsigned int jj = jj_lower; jj <= jj_upper; jj++)
                {
                    if (fast){
                        pos(0) = binners[0].get_bin_midpoint(ii);
                        pos(1) = binners[1].get_bin_midpoint(jj);
                        double corr_0 =  (calc_boundary_correction_term(pos(0), corr_stddev_0, binners[0].get_min_val(), binners[0].get_max_val() ));
                        double corr_1 =  (calc_boundary_correction_term(pos(1), corr_stddev_1, binners[1].get_min_val(), binners[1].get_max_val() ));
                        double val =  (this->bin_volume * (( pdf.getPdf(pos) ))/(corr_0*corr_1));
                        add_count(make_indices(ii,jj), norm_term * weight*(val));
                    }
                    else {
                        pos_UU(0) = binners[0].get_bin_upperbound(ii);
                        pos_UU(1) = binners[1].get_bin_upperbound(jj);
                        pos_UL(0) = binners[0].get_bin_upperbound(ii);
                        pos_UL(1) = binners[1].get_bin_lowerbound(jj);
                        pos_LU(0) = binners[0].get_bin_lowerbound(ii);
                        pos_LU(1) = binners[1].get_bin_upperbound(jj);
                        pos_LL(0) = binners[0].get_bin_lowerbound(ii);
                        pos_LL(1) = binners[1].get_bin_lowerbound(jj);
                        double corr_0 =(calc_boundary_correction_term(pos_UU(0), corr_stddev_0, binners[0].get_min_val(), binners[0].get_max_val() )
                                    + calc_boundary_correction_term(pos_UL(0), corr_stddev_0, binners[0].get_min_val(), binners[0].get_max_val() )
                                    + calc_boundary_correction_term(pos_LU(0), corr_stddev_0, binners[0].get_min_val(), binners[0].get_max_val() )
                                    + calc_boundary_correction_term(pos_LL(0), corr_stddev_0, binners[0].get_min_val(), binners[0].get_max_val() )) / 4.0;

                        double corr_1 =(calc_boundary_correction_term(pos_UU(1), corr_stddev_0, binners[1].get_min_val(), binners[1].get_max_val() )
                                    + calc_boundary_correction_term(pos_UL(1), corr_stddev_0, binners[1].get_min_val(), binners[1].get_max_val() )
                                    + calc_boundary_correction_term(pos_LU(1), corr_stddev_0, binners[1].get_min_val(), binners[1].get_max_val() )
                                    + calc_boundary_correction_term(pos_LL(1), corr_stddev_0, binners[1].get_min_val(), binners[1].get_max_val() )) / 4.0;
                        double val =  (this->bin_volume * ((pdf.getPdf(pos_UU) + pdf.getPdf(pos_UL) + pdf.getPdf(pos_LU) + pdf.getPdf(pos_LL))/4.0 ))/(corr_0*corr_1);
                        add_count(make_indices(ii,jj), norm_term * weight*(val));
                    }
                    //total_add+= weight*(val);
                }
            }
        }


    }

	void add_gaussian_count (const double weight, double_array_type p, double_array_type stddev, const double stddev_cutoff = 3.0, const bool fast = true){
	    new_add_gaussian_count(weight, p , stddev, stddev_cutoff, fast);
	    return;


        /*
        uint_array_type indicies = get_indicies(p);
        if (stddev[0] == 0 && stddev[1] == 0){
            //std::cout << "adding count " << weight << " to: " << indicies[0] << " " << indicies[1] << " current_count: " << get_count(indicies) << std::endl;
            add_count(indicies, weight);
            //std::cout << "new count: " << get_count(indicies) << std::endl;
            return;
        }

        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar = Eigen::MatrixXd::Zero(2, 2);
        Eigen::Matrix<double,Eigen::Dynamic, 1> mean = Eigen::MatrixXd::Zero(2, 1);

        covar(0,0) = stddev[0] * stddev[0];
        covar(1,1) = stddev[1] * stddev[1];
        mean(0) = p[0];
        mean(1) = p[1];

        MultivariateNormalPDF pdf;
        pdf.setMean(mean);
        pdf.setCovar(covar);

        const double cutoff_range_0 = stddev_cutoff * stddev[0];
        const double cutoff_range_1 = stddev_cutoff * stddev[1];

        const unsigned int cutoff_bin_range_0 = std::ceil(cutoff_range_0 / binners[0].get_bin_size());
        const unsigned int cutoff_bin_range_1 = std::ceil(cutoff_range_1 / binners[1].get_bin_size());

        const unsigned int ii_lower = cutoff_bin_range_0 <= indicies[0] ? indicies[0] - cutoff_bin_range_0 : 0;
        const unsigned int jj_lower = cutoff_bin_range_1 <= indicies[1] ? indicies[1] - cutoff_bin_range_1 : 0;

        const unsigned int ii_upper = (indicies[0] + cutoff_bin_range_0 + 1 ) < num_bins[0] ? indicies[0] + cutoff_bin_range_0 + 1 : num_bins[0];
        const unsigned int jj_upper = (indicies[1] + cutoff_bin_range_1 + 1 ) < num_bins[1] ? indicies[1] + cutoff_bin_range_1 + 1 : num_bins[1];

        //const double val_at_mean = pdf.getPdf(mean);

        //std::cout << "ii range: " << ii_lower << "\t" << ii_upper << std::endl;
        //std::cout << "jj range: " << jj_lower << "\t" << jj_upper << std::endl;


        Eigen::Matrix<double,Eigen::Dynamic, 1> pos = Eigen::MatrixXd::Zero(2, 1);

        //std::cout << "bin_vol =  " << bin_volume << std::endl;


        const double norm_term = 1;//weight / total_add;

        for (unsigned int ii = ii_lower; ii < ii_upper; ii++){
            for (unsigned int jj = jj_lower; jj < jj_upper; jj++){
                pos(0) = binners[0].get_bin_midpoint(ii);
                pos(1) = binners[1].get_bin_midpoint(jj);
                double val = norm_term * this->bin_volume * pdf.getPdf(pos);
                add_count(make_indices(ii,jj), weight*(val));
                //total_add+= weight*(val);

            }
        }

        //std::cout << "total_add: " << total_add << std::endl;
        */

	}


	virtual std::ostream& output(std::ostream& output) const{


	    //double sum_vals = 0;



        /*
	    for (multi_array_type::const_iterator it = hist.begin(); it != hist.end(); it++){
            sum_vals += *it;
            sum_vals += std::exp(-(*it);
	    }
	    */

        //auto ntup = gen_tuple<std::vector<unsigned long>, 2 >();
        std::vector<unsigned int> this_vec0 = std::vector<unsigned int>(num_bins[0], 0);
        std::iota(this_vec0.begin(), this_vec0.end(), 0);

        std::vector<unsigned int> this_vec1 = std::vector<unsigned int>(num_bins[1], 0);
        std::iota(this_vec1.begin(), this_vec1.end(), 0);

        //auto ntup = std::make_tuple(this_vec0, this_vec1);


        double max_value = std::numeric_limits< double>::lowest();
        for (auto&& t : iter::product(this_vec0, this_vec1)) {
            uint_array_type pos = {std::get<0>(t), std::get<1>(t)};
            const double val_here = get_count(pos);
            if (val_here > max_value){
                max_value = val_here;
            }
        }

        const double avg_val = running_sum / double(hist.num_elements());
        double sum_exp_vals = 0;
        for (auto&& t : iter::product(this_vec0, this_vec1)) {
            uint_array_type pos = {std::get<0>(t), std::get<1>(t)};
            //sum_vals += get_count(pos);
            sum_exp_vals += std::exp(get_count(pos) - max_value);
        }

        //const double avg_val = sum_vals / double(hist.num_elements());

        //std::cout << "# DEBUG" << this_vec0.size() << std::endl;
        for (auto&& t : iter::product(this_vec0, this_vec1)) {
            //std::cout << "# BLAH" << std::endl;
            uint_array_type pos = {std::get<0>(t), std::get<1>(t)};
            output << std::get<0>(t) << "\t"
                    << std::get<1>(t)<< "\t"
                    << binners[0].get_bin_midpoint(std::get<0>(t)) << "\t"
                    << binners[1].get_bin_midpoint(std::get<1>(t)) << "\t"
                    << get_count(pos) << "\t"
                    << get_count(pos) - avg_val << "\t"
                    << std::exp(get_count(pos) - max_value) / sum_exp_vals << "\t"
                    << get_count(pos) / running_sum << "\t"
                    << std::endl;
        }


        return output;
	}

    virtual std::ostream& output(std::ostream& output, const double minus_factor, const double sum_exp_vals, const bool ignore_bin_volume = false) const{
                std::vector<unsigned int> this_vec0 = std::vector<unsigned int>(num_bins[0], 0);
        std::iota(this_vec0.begin(), this_vec0.end(), 0);

        std::vector<unsigned int> this_vec1 = std::vector<unsigned int>(num_bins[1], 0);
        std::iota(this_vec1.begin(), this_vec1.end(), 0);

        //auto ntup = std::make_tuple(this_vec0, this_vec1);



        /*
        double max_value = std::numeric_limits< double>::lowest();
        for (auto&& t : iter::product(this_vec0, this_vec1)) {
            uint_array_type pos = {std::get<0>(t), std::get<1>(t)};
            const double val_here = get_count(pos);
            if (val_here > max_value){
                max_value = val_here;
            }
        }

        const double avg_val = running_sum / double(hist.num_elements());
        double sum_exp_vals = 0;
        for (auto&& t : iter::product(this_vec0, this_vec1)) {
            uint_array_type pos = {std::get<0>(t), std::get<1>(t)};
            //sum_vals += get_count(pos);
            sum_exp_vals += std::exp(get_count(pos) - max_value);
        }
        */

        //const double avg_val = sum_vals / double(hist.num_elements());

        //std::cout << "# DEBUG" << this_vec0.size() << std::endl;
        for (auto&& t : iter::product(this_vec0, this_vec1)) {
            //std::cout << "# BLAH" << std::endl;
            uint_array_type pos = {std::get<0>(t), std::get<1>(t)};
            output << std::get<0>(t) << "\t"
                    << std::get<1>(t)<< "\t"
                    << binners[0].get_bin_midpoint(std::get<0>(t)) << "\t"
                    << binners[1].get_bin_midpoint(std::get<1>(t)) << "\t"
                    << get_count(pos) << "\t"
                    << get_count(pos) + (ignore_bin_volume ? 0 : std::log(this->bin_volume)) - minus_factor << "\t"
                    << std::exp(get_count(pos) + (ignore_bin_volume ? 0 : std::log(this->bin_volume)) - minus_factor) / sum_exp_vals << "\t"
                    << get_count(pos) / running_sum << "\t"
                    << std::endl;
        }


        return output;
	}


};



#endif
