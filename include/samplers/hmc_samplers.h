// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_SAMPLERS_H
#define HMC_SAMPLERS_H

#include <cmath>

template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier(Polytope &P, Point &p, PointList randPoints, NT &a, int n, int N) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);

    MT A = P.get_mat();
    VT b = P.get_vec();
    unsigned int d = P.dimension();
    unsigned int m = A.nrows();

    VT v0(d), x0(d);
    for (int i = 0; i < d; ++i) {
        v0(i) = rdist(rng);
        x0(i) = p[i];
    }

    s0 = A*x0;
    sv0 = A*v0;
    MT M = - A*A.transpose();

    VT c(3);
    cj(0) = 0.0; cj(1) = 0.587785; cj(2) = 0.951056; //Chebyshev nodes

    MT AA = MT::Ones(n+1, n+1);
    MT T = MT::Zero(n+1,n+1);
    VT S = VT::Zero(n+1);

    for (int j = 1; j < n+1; ++j) {
        AA.col(j) = cj.pow(j);
        S(j) = NT(j)*cj,pow(j-1);
        if (j>1) T.col(j) = NT(j)*NT(j-1)*cj.pow(j-2);
    }
    AAinv = A.inverse();
    T = T*AAinv;
    S = S*AAinv;
    MT pinvA = A.completeOrthogonalDecomposition().pseudoInverse()

    for (int i = 0; i < N; ++i) {

        get_next_hmc_logbarrier<RNGType>(x0, T, S, M, s0, sv0, cj, pinvA, b, a);

        for (int i = 0; i < d; ++i) {
            v0(i) = rdist(rng);
            p[i] = x0(i);
        }
        randPoints.push_back(p);
        s0 = A*x0;
        sv0 = A*v0;
    }

}


template <class RNGType, class VT, class MT, typename NT>
void get_next_hmc_logbarrier(VT &x0, MT &T, VT &S, MT &M, VT &s0, VT &sv0, VT &cj, MT pinvA, VT &b, NT &a) {

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);

    unsigned int m = s0.size();
    unsigned int d = x0.size();
    unsigned int n = cj.size()-1;

    MT J = MT::Zero(m*(n+1), m*(n+1));
    VT F = VT::Zero(M*(n+1));
    VT s1(m*(n+1));
    for (int j = 0; j < m*(n+1); ++j) {
        s1 = rdist(rng);
    }
    VT s2 = 10*s1;

    for (int i = 0; i < m; ++i) {
        J(m * (n - 1) + i, i * (n + 1)) = 1.0;
        for (int j = 0; j < n + 1; ++j) {
            J(m * n + i, (i - 1) * (n + 1) + j) = S(1, j);
        }

        for (int j = 1; j < n - 1; ++j) {
            for (int k = 0; k < n + 1; ++k) {
                J((i - 1) * (n - 1) + j - 1, (i - 1) * (n + 1) + k) = T(j, k);
            }
        }
    }

    while ((s1-s2).abs().MaxCoeff()>=0.000001) {

        for (int i = 0; i < m; ++i) {

            F(m * (n - 1) + i) = s1((i - 1) * (n + 1)) - s0(i);
            F(m * n + i) = S.dot(s1.segment(((i - 1) * (n + 1)), ((i - 1) * (n + 1) + n))) - sv0(i);

            for (int j = 1; j < n - 1; ++j) {
                J((i - 1) * (n - 1) + j - 1, (i - 1) * (n + 1) + j) = T(j, j) -
                        M(i, i) * a * (1 / ((b(i) - s1((i - 1) * (n + 1) + j)) * (b(i) - s1((i - 1) * (n + 1) + j))));
                for (int k = 0; k < m; ++k) {
                    vec(k) = a / (b(k) - s1((k - 1) * (n + 1) + j));
                }
                F((i - 1) * (n - 1) + j - 1) =
                        T.row(j).dot(s1.segment(((i - 1) * (n + 1)), ((i - 1) * (n + 1) + n))) - M.row(i).dot(vec);
            }
        }

        update_svals(s1, s2, J, F);

    }
    

}


#endif
