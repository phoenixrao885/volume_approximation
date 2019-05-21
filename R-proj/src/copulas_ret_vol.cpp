// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "simplex_samplers2.h"
#include "families.h"
#include "copulas_hnr.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix get_copulas (Rcpp::NumericMatrix RetMat, Rcpp::NumericMatrix EllMats,
                             Rcpp::Nullable<unsigned int> Win = R_NilValue,
                             Rcpp::Nullable<unsigned int> M = R_NilValue,
                             Rcpp::Nullable<unsigned int> N = R_NilValue,
                             Rcpp::Nullable<double> error = R_NilValue,
                             Rcpp::Nullable<double> prob = R_NilValue){

    typedef double NT;
    typedef boost::mt19937 RNGType;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    unsigned int MM = (M.isNotNull()) ? Rcpp::as<unsigned int>(M) : 100;
    unsigned int NN = (N.isNotNull()) ? Rcpp::as<unsigned int>(N) : 500000;
    unsigned int W = (Win.isNotNull()) ? Rcpp::as<unsigned int>(Win) : 60;
    NT e = (error.isNotNull()) ? Rcpp::as<double>(error) : 0.1;
    NT pr = (prob.isNotNull()) ? Rcpp::as<double>(prob) : 0.9;

    int d = EllMats.cols();
    int K = EllMats.rows() / d;
    std::cout<<EllMats.cols()<<" "<<EllMats.rows()<<"\n"<<std::endl;
    std::cout<<"K = "<<K<<" d = "<<d<<std::endl;

    std::vector<MT> ellipsoids;
    MT hyperplanes(K, d);
    get_bodies<VT, NT>(Rcpp::as<MT>(RetMat), Rcpp::as<MT>(EllMats), ellipsoids,  hyperplanes, W);
    std::cout<<"bodies ok"<<std::endl;
    std::cout<<hyperplanes.cols()<<" "<<hyperplanes.rows()<<"\n"<<std::endl;
    //for (int k = 0; k < 2; ++k) {
        std::cout<<ellipsoids[0].cols()<<" "<<ellipsoids[0].rows()<<"\n"<<std::endl;
    //}
    std::vector<MT> copulas = get_copulas_uniform<RNGType, VT> (ellipsoids, hyperplanes, e, pr, MM, NN);

    MT ret_copulas(K*MM,MM);
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < MM; ++j) {
            ret_copulas.row(j+i*MM) = copulas[i].row(j);
        }
    }

    return Rcpp::wrap(ret_copulas);
}