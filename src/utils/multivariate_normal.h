#ifndef MULTIVARIATE_NORMAL_H_INCLUDED
#define MULTIVARIATE_NORMAL_H_INCLUDED

#ifndef M_PI
#define M_PI REAL(3.1415926535897932384626433832795029)
#endif

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
#include <vector>
#include  <cmath>

#define MAXBUFSIZE  ((int) 1e6)




class MultivariateNormalPDF{
private:

    Eigen::Matrix<double,Eigen::Dynamic,1> mean ;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar_inv;

    double det, logdet;
    double normalising_term, log_normalising_term;

public:




    MultivariateNormalPDF(){
    }

    Eigen::Matrix<double,Eigen::Dynamic,1> get_mean() const{
        return mean;
    }

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> get_covar() const{
        return covar;
    }

    void setCovar(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar_, const bool use_factorization = true){
        covar = covar_;

        if (use_factorization == false){
            covar_inv = covar.inverse();
            det = covar.determinant();
            normalising_term = 1.0 / std::sqrt(std::pow(2.0 * M_PI, covar.rows()) * det);
            log_normalising_term = std::log(normalising_term);
        }
        else {

            Eigen::FullPivHouseholderQR< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > qr(covar);
            covar_inv = qr.inverse();
            logdet = qr.logAbsDeterminant();
            det = std::exp(logdet);
            log_normalising_term = 0.5 * (-std::log(std::pow(2.0 * M_PI, covar.rows())) - logdet);
            normalising_term = std::exp(log_normalising_term);
        }
    }

    void setMean(Eigen::Matrix<double,Eigen::Dynamic,1> mean_){
        mean = mean_;
    }

    double  getPdf(Eigen::Matrix<double,Eigen::Dynamic,1> val){
        /*
        Eigen::Matrix<double,Eigen::Dynamic,1>  val_cent = val - mean;
        const double exp_term =  std::exp((-0.5 * (val_cent.transpose()* covar_inv * val_cent))(0,0));
        const double prob = normalising_term * exp_term;
        */
        return std::exp(this->getLogPdf(val)); //prob;
    }


    double  getLogPdf(Eigen::Matrix<double,Eigen::Dynamic,1> val){
        Eigen::Matrix<double,Eigen::Dynamic,1>  val_cent = val - mean;
        //std::cout << val_cent << std::endl;
        const double exp_term =  (-0.5 * (val_cent.transpose()* covar_inv * val_cent))(0,0);
        //std::cout << "exp_term: " << exp_term << std::endl;
        const double log_prob = log_normalising_term + exp_term;
        //std::cout << "log_normalising_term " << log_normalising_term << std::endl;
        //std::cout << "normalising_term " << normalising_term << std::endl;
        //std::cout << "det " << det << std::endl;
        return log_prob;
    }

    double  getLogPdfNoNorm(Eigen::Matrix<double,Eigen::Dynamic,1> val){
        Eigen::Matrix<double,Eigen::Dynamic,1>  val_cent = val - mean;
        //std::cout << val_cent << std::endl;
        const double exp_term =  (-0.5 * (val_cent.transpose()* covar_inv * val_cent))(0,0);
        //std::cout << "exp_term: " << exp_term << std::endl;
        const double log_prob = exp_term;
        //std::cout << "log_normalising_term " << log_normalising_term << std::endl;
        //std::cout << "normalising_term " << normalising_term << std::endl;
        //std::cout << "det " << det << std::endl;
        return log_prob;
    }

};


#endif // MULTIVARIATE_NORMAL_H_INCLUDED
