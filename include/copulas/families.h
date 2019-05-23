 // VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

#ifndef FAMILIES_H
#define FAMILIES_H

#include "exact_vols.h"


 template <class VT, typename NT, class MT>
 MT get_hyp_consts(MT &hyperplanes, int M) {

     int d = hyperplanes.cols(), K = hyperplanes.rows(), con;
     VT hyp(d);
     NT cmin, cmax, min, max, med, ratio, error = 0.001, curr_ratio = 1.0 / NT(M);
     MT constants(K, M - 1);

     for (int i = 0; i < K; ++i) {

         hyp = hyperplanes.row(i);
         //std::cout << hyp << std::endl;
         cmin = hyp.minCoeff();
         cmax = hyp.maxCoeff();
         min = cmin;
         max = cmax;

         for (int j = 0; j < M - 1; ++j) {

             while (true) {

                 med = 0.5 * (min + max);
                 ratio = vol_Ali(hyp, -med);
                 if (ratio - j * curr_ratio > curr_ratio &&
                     (ratio - j * curr_ratio - curr_ratio) / curr_ratio > error) {
                     max = med;
                 } else if (ratio - j * curr_ratio < curr_ratio &&
                            (curr_ratio - (ratio - j * curr_ratio)) / curr_ratio > error) {
                     min = med;
                 } else {
                     constants(i, j) = med;
                     break;
                 }
             }

             min = med;
             max = cmax;

         }

     }
     return constants;
 }


template <class VT, typename NT, class MT>
void get_bodies (MT RetMat, MT allEll, std::vector<MT> &ellipsoids,  MT &hyperplanes, int W) {

    int d = RetMat.cols(), L = RetMat.rows();
    VT hyp(d);
    MT compmat(W, d);
    MT temp_ell(d,d);

    for (int i = 0; i < L-W+1; ++i) {

        hyp = VT::Ones(d);
        for (int k = 0; k < d; ++k) {
            for (int j = 0; j < W; ++j) {
                hyp(k) = hyp(k) * (1.0 + RetMat(i+j, k));
            }
            hyp(k) = hyp(k) - 1.0;
        }
        hyperplanes.row(i) = hyp;

        for (int l = 0; l < d; ++l) {
            temp_ell.row(l) = allEll.row(i*d+l);
        }
        ellipsoids.push_back(temp_ell);
    }
}



template <class RNGType, class VT, class MT, typename NT>
std::pair<MT, MT> get_constants(std::vector<MT> &ellipsoids,  MT &hyperplanes, std::vector<MT> &copulas, int M, int N,
                                std::vector<NT> &mins) {

    MT copula(M, M);
    int d = hyperplanes.row(0).size(), K = ellipsoids.size(), row, col;
    MT points(d, N);
    exp_simplex<NT, RNGType>(d, N, points);
    //std::cout<<points.cols()<<" "<<points.rows()<<"\n"<<std::endl;


    MT hyp_vals = hyperplanes * points;
    MT hyp_cons = get_hyp_consts<VT, NT>(hyperplanes, M);
    MT ell_cons(K, M-1), vecs(d, N), ell_vals(K, N);
    VT row_vals(N);

    int count_ell = 0;
    NT pos = 1.0 / NT(M);
    for (typename std::vector<MT>::iterator ellit = ellipsoids.begin();  ellit!=ellipsoids.end(); ++ellit, ++count_ell) {

        vecs = (*ellit) * points;
        for (int i = 0; i < N; ++i) ell_vals(count_ell, i) = points.col(i).dot(vecs.col(i));

        row_vals = ell_vals.row(count_ell);
        std::sort(row_vals.data(), row_vals.data() + row_vals.size());

        for (int i=1; i<M; i++) {
            ell_cons(count_ell, i - 1) = row_vals(((int) std::floor(i * (pos) * (NT(N)))));
        }

        copula = MT::Zero(M,M);
        for (int k = 0; k < N; ++k) {

            row = -1;
            col = -1;
            for (int j = 0; j < M - 1; ++j) {
                if (ell_vals(count_ell, k) < ell_cons(count_ell, j)) {
                    row = j;
                    break;
                }
            }
            for (int j = 0; j < M - 1; ++j) {
                if (hyp_vals(count_ell, k) < hyp_cons(count_ell, j)) {
                    col = j;
                    break;
                }
            }
            if (col == -1) col = M - 1;
            if (row == -1) row = M - 1;
            copula(row, col) = copula(row, col) + 1.0;
        }
        mins.push_back(copula.minCoeff()/NT(N));
        copulas.push_back(copula);
        std::cout<<"number of copulas = "<<copulas.size()<<std::endl;
    }

    return std::pair<MT, MT> (ell_cons, hyp_cons);
}


#endif
