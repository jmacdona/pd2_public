#ifndef HIST_UTILS_H_INCLUDED
#define HIST_UTILS_H_INCLUDED

#include <cmath>
#include <vector>
#include "boost/multi_array.hpp"
#include <iostream>
#include <array>
#include <numeric>
#include <utility>
#include "product.hpp"


//! This class determines the corresponding bin index for a given value in one dimension
class hist_binner{

public:
        unsigned int num_bins;
        double min_val, max_val, bin_size;

public:
    hist_binner(unsigned int num_bins_, const double min_val_, const double max_val_) : num_bins(num_bins_),
        min_val(min_val_), max_val(max_val_), bin_size( (max_val_ - min_val_)/double(num_bins_) ) {

    }

    hist_binner(){
    }

    //! this returns the bin number where low_bin_boundary <= val < upper_bin_boundary except for the last bin where low_bin_boundary <= val in order to guarantee returned bins fall in a defined range
    unsigned int get_bin(double val) const{
        unsigned int bin = (unsigned int)std::floor( (val - min_val) / bin_size );
        return bin < num_bins ? bin  : (num_bins-1);
    }

	double get_bin_lowerbound(unsigned int bin) const{
		return (double(bin) * bin_size) + min_val;
	}

	double get_bin_midpoint(unsigned int bin) const{
		return (double(bin) * bin_size) + (bin_size/2.0)  + min_val;
	}

	double get_bin_upperbound(unsigned int bin) const{
		return (double(bin+1) * bin_size) + min_val;
	}

	double get_bin_size() const{
		return bin_size;
	}
	double get_min_val() const {
		return min_val;
	}
	double get_max_val() const {
		return max_val;
	}

	unsigned int get_num_bins() const{
		return num_bins;
	}


};


/*
template<size_t...>
struct index_sequence{};
namespace details {
  template<size_t count, size_t...Is>
  struct mis_helper:mis_helper<count-1, count-1, Is...> {};
  template<size_t...Is>
  struct mis_helper<0,Is...> {
    using type=index_sequence<Is...>;
  };
}
template<size_t count>
using make_index_sequence=typename details::mis_helper<count>::type;

template<typename Func, typename Tup, std::size_t... index>
decltype(auto) invoke_helper(Func&& func, Tup&& tup, index_sequence<index...>)
{
    return func(std::get<index>(std::forward<Tup>(tup))...);
}

template<typename Func, typename Tup>
decltype(auto) invoke(Func&& func, Tup&& tup)
{
    constexpr auto Size = std::tuple_size<typename std::decay<Tup>::type>::value;
    return invoke_helper(std::forward<Func>(func),
                         std::forward<Tup>(tup),
                         make_index_sequence<Size>{});
}


template<std::size_t> struct int_{};

template <class Tuple, size_t Pos>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<Pos> ) {
  out << std::get< std::tuple_size<Tuple>::value-Pos >(t) << '\t';
  return print_tuple(out, t, int_<Pos-1>());
}

template <class Tuple>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<1> ) {
  return out << std::get<std::tuple_size<Tuple>::value-1>(t);
}

template <class... Args>
std::ostream& operator<<(std::ostream& out, const std::tuple<Args...>& t) {
  print_tuple(out, t, int_<sizeof...(Args)>());
}
*/


/*
template<size_t, class T>
using T_ = T;

template<class T, size_t... Is>
decltype(auto) gen_tuple(index_sequence<Is...>) { return std::tuple<T_<Is, T>...>{}; }

template<class T, size_t N>
decltype(auto) gen_tuple() { return gen_tuple<T>(make_index_sequence<N>{}); }
*/


template <unsigned int N>
class multi_hist{
public:





    typedef boost::multi_array<double, N> multi_array_type;
    typedef boost::array< unsigned int, N > uint_array_type;
    typedef boost::array< int, N > int_array_type;
    typedef boost::array< hist_binner, N > hist_binner_array_type;
    typedef boost::array< double, N > double_array_type;


    //! the number of bins in each dimension
    uint_array_type num_bins;
    //! array of binners for each dimension
    hist_binner_array_type binners;
    //! minimum and maximum values for each dimension
    double_array_type min_vals, max_vals;

    //! stores the actual counts/values at each multidimensional bin location
    multi_array_type hist;

    double incr;
    double running_sum;
    double bin_volume;


public:


    // TODO: make this variadic
    static uint_array_type make_indices(const unsigned int i, const unsigned int j){
    	uint_array_type indices = {i, j};
    	return indices;
    }

    multi_hist() : num_bins(), binners(), min_vals(), max_vals(), hist(), incr(0), running_sum(0), bin_volume(0){
    }

    multi_hist(const multi_hist& orig_) : num_bins(orig_.num_bins), binners(orig_.binners), min_vals(orig_.min_vals), max_vals(orig_.max_vals), incr(orig_.incr), running_sum(orig_.running_sum), bin_volume(orig_.bin_volume)  {
        // copy hist
        hist.resize(num_bins);
        hist = orig_.hist;
    }

    multi_hist& operator=(const multi_hist& orig_) {

        num_bins = orig_.num_bins;
        binners = orig_.binners;
        min_vals = orig_.min_vals;
        max_vals = orig_.max_vals;
        incr = orig_.incr;
        running_sum = orig_.running_sum;
        bin_volume =orig_.bin_volume;

        hist.resize(num_bins);
        hist = orig_.hist;

        return *this;
    }


    multi_hist(uint_array_type num_bins_, double_array_type min_vals_, double_array_type max_vals_, double incr_, double pseudo_count ) :
        num_bins(num_bins_), min_vals(min_vals_), max_vals(max_vals_), hist(num_bins), incr(incr_), running_sum(pseudo_count * hist.num_elements())   {


        //std::cout << "#\n#initialised hist with dimensionality" << std::endl;
        //std::cout << "#" << num_bins_.size() << std::endl;
        //std::cout << "#" << hist.dimensionality << std::endl;
        //std::cout << "#" << hist.num_elements() << std::endl;
        //hist = multi_array_type(num_bins);
        //std::cout << "#filling hist with zeros" << std::endl;
        std::fill(hist.data(), hist.data() + hist.num_elements(), pseudo_count);
        //std::cout << "#init hist binners" << std::endl;
        bin_volume = 1;
        for (unsigned int ii = 0; ii < N; ii++){
            binners[ii] = hist_binner(num_bins[ii], min_vals[ii], max_vals[ii]);
            bin_volume *= binners[ii].get_bin_size();
        }

        //std::cout << "#done." << std::endl;
    }

	uint_array_type get_indicies(double_array_type p) const{
		uint_array_type indicies;
		for (unsigned int ii = 0; ii < N; ii++){
			indicies[ii] = binners[ii].get_bin(p[ii]);
		}

		return indicies;
	}


    void add_count(uint_array_type indicies, const double val){
        running_sum += val;
        hist(indicies) += val;
    }

    void add_count(double_array_type p, const double val){
        add_count(get_indicies(p), val);
        //hist(get_indicies(p)) += val;
    }

    // add counts to histogram
    void add_count(double_array_type p){
        add_count(p, incr);
    }

	double get_count(double_array_type p) const{
		return hist(get_indicies(p));
	}

    uint_array_type get_num_bins() const{
        return num_bins;
    }

    double get_count(uint_array_type pos) const{
        return hist(pos);
    }

	void set_count(uint_array_type pos, const double val) {
	    running_sum -= hist(pos);
	    running_sum += val;
		hist(pos) = val;
	}

	void set_all_counts(double val){
        std::fill(hist.data(), hist.data() + hist.num_elements(), val);
        running_sum = val * hist.num_elements();
	}

	double get_running_sum() const{
	    return running_sum;
	}

	double get_vol_scaled_max_value(const bool ignore_volume = false) const {
	    double max_value = std::numeric_limits< double>::lowest();
	    for (unsigned int ii = 0; ii < hist.num_elements(); ii++){
                if (*(hist.data()+ii) > max_value){
                    max_value = *(hist.data()+ii)+ (ignore_volume ? 0 : std::log(this->bin_volume));
                }
	    }
	    return max_value;
	}

	double get_vol_scaled_sum_exp(const double minus_factor, const bool ignore_volume = false) const {
	    double sum = 0;
	    for (unsigned int ii = 0; ii < hist.num_elements(); ii++){
            sum += std::exp((*(hist.data()+ii) - minus_factor + (ignore_volume ? 0 : std::log(this->bin_volume))));
	    }
        return sum;
	}

	double get_sum() const {
	    double sum = 0;
	    for (unsigned int ii = 0; ii < hist.num_elements(); ii++){
            sum += *(hist.data()+ii);
	    }
        return sum;
	}

	double get_bin_volume() const{
	    return this->bin_volume;
	}


	virtual std::ostream& output(std::ostream& output) const{


        output << "hist_utils.output(): ERROR: not currently supported" << std::endl;

        /*
	    auto ntup = gen_tuple<std::vector<unsigned long>, N >();

	    std::array< std::vector<unsigned long>, N > narray;
	    for (int dim = 0 ; dim < hist.dimensionality; dim++){
            std::vector<unsigned long> this_vec = std::vector<unsigned long>(0, num_bins[dim]);
            std::iota(this_vec.begin(), this_vec.end(), 0);
            narray[dim] = this_vec;
	    }
	    */


        /*
        // this is probably not going to work as need to unpack narray
	    for (auto&& t : invoke(iter::product, ntup)) {
            output << t << std::endl;
        }
        */

        return output;
	}

};


#endif // HIST_UTILS_H_INCLUDED
