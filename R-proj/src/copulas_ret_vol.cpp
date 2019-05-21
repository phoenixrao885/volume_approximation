// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix copulas (Rcpp::NumericMatrix RetMat, Rcpp::NumericMatrix EllMats,
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
    unsigned int NN = (N.isNotNull()) ? Rcpp::as<unsigned int>(N) : 100000;
    unsigned int W = (Win.isNotNull()) ? Rcpp::as<unsigned int>(Win) : 60;
    NT e = (error.isNotNull()) ? Rcpp::as<double>(error) : 0.1;
    NT pr = (prob.isNotNull()) ? Rcpp::as<double>(prob) : 0.9;

    int d = EllMats.ncol();
    int K = EllMats.nrow() / d;

    std::vector<MT> ellipsoids;
    MT hyperplanes(K, d);
    get_bodies (Rcpp::as<MT>(RetMat), Rcpp::as<MT>(EllMats), ellipsoids,  hyperplanes, W);
    std::vector<MT> copulas = copulas_uniform (ellipsoids, hyperplanes, e, pr, MM, NN);

    MT ret_copulas(K*MM,MM);
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < MM; ++j) {
            ret_copulas.row(j+i*MM) = copulas[i].row(j);
        }
    }

    return Rcpp::wrap(ret_copulas);
}