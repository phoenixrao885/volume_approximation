// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef HZONO_VOL_H
#define HZONO_VOL_H

#include "ball_annealingGl.h"
#include "hpoly_annealing.h"
#include "esti_ratioGl.h"
#include "ZonoIntersectHPoly.h"


template <class Hpolytope, class Zonotope, class UParameters, class AParameters, class GParameters, class Point, typename NT>
NT vol_hzono (Zonotope &ZP, UParameters &var, AParameters &var_ban, GParameters &var_g,
              std::pair<Point,NT> &InnerB, NT &nballs, Rcpp::Function diam_zono) {

    typedef typename Zonotope::VT VT;
    typedef typename Zonotope::MT MT;
    typedef typename UParameters::RNGType RNGType;

    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, rmax = var_ban.rmax,
            radius = InnerB.second, round_value = 1.0, e = var.error, ratio, vol, alpha = var_ban.alpha;
    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;

    //var.verbose = true;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;

    MT V = ZP.get_mat();
    MT G = V.transpose();
    int m = G.cols();
    std::list<Point> randPoints;

    MT AA; VT b;
    AA.resize(2 * m, m);
    b.resize(2 * m);
    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < m; ++j) {
            if (i == j) {
                AA(i, j) = 1.0;
            } else {
                AA(i, j) = 0.0;
            }
        }
    }
    for (unsigned int i = 0; i < m; ++i) {
        b(i + m) = 1.0;
        for (unsigned int j = 0; j < m; ++j) {
            if (i == j) {
                AA(i + m, j) = -1.0;
            } else {
                AA(i + m, j) = 0.0;
            }
        }
    }

    //MT A3;
    //NT volh;
    MT T = ZP.get_T();
    MT Tt = T.transpose();
    MT A2 = AA * Tt;
    MT B = G * Tt;
    MT A3 = A2 * B.inverse();

    MT A = A3*G;

    //std::cout<<"M*G = "<<A3*G<<std::endl;

    Hpolytope HP;
    HP.init(n,A3,b);

    VT Zs_max(2*m);
    UParameters var3 = var;
    var3.cdhr_walk = true;
    var3.ball_walk = false;
    var3.rdhr_walk = false;
    var3.bill_walk = false;
    get_hdelta(ZP, HP, Zs_max, lb, ub, ratio, var3);
    Hpolytope HP2=HP;

    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();

    typedef ZonoIntersectHPoly<Zonotope , Hpolytope > ZonoHP;
    std::vector<Hpolytope > HPolySet;
    std::vector<NT> ratios;
    NT p_value = 0.1;
    ZonoHP zb1, zb2;

    //NT ole = 0.6*Rcpp::as<NT>(diam_zono(Rcpp::wrap(ZP.get_mat().transpose()), Rcpp::wrap(A), Rcpp::wrap(HP.get_vec())));
    //std::cout<<"ole = "<<ole<<std::endl;

    get_sequence_of_zonopolys<ZonoHP>(ZP, HP2, HPolySet, Zs_max, ratios, N*nu, nu, lb, ub, alpha, var);
    nballs = NT(HPolySet.size()+1);

    int mm=HPolySet.size()+2;
    int mm2=mm+1;
    prob = std::pow(prob, 1.0/NT(mm2));
    NT er0 = e/(2.0*std::sqrt(NT(mm2)));
    NT er1 = (e*std::sqrt(2.0*NT(mm2)-1))/(std::sqrt(2.0*NT(mm2)));
    NT Her = e/(2.0*std::sqrt(NT(mm2)));

    var_g.error = Her/2.0;
    NT nballs2;
    vol = volume_gaussian_annealing(HP, var_g, var, InnerBall,nballs2);

    if(verbose) std::cout<<"\nvol of h-polytope = "<<vol<<"\n"<<std::endl;
    if (!window2) {
        UParameters var2 = var;
        var2.cdhr_walk = true;
        var2.ball_walk = false;
        var2.rdhr_walk = false;
        var2.bill_walk = false;
        var2.walk_steps = 10+2*n;
        vol *= esti_ratio_interval<RNGType, Point>(HP2, ZP, ratio, er0, win_len, N*nu, prob, var2);
    } else {
        vol *= esti_ratio<RNGType, Point>(HP2, ZP, ratio, er0, var_g.W, N*nu, var);
    }

    Hpolytope b1, b2;
    if (HPolySet.size()==0) {
        if (verbose) std::cout << "no hpoly | ratio = " << ratios[0] << std::endl;
        if (ratios[0]!=1) {
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(ZP, HP2, ratios[0], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(ZP, HP2, ratios[0], er1, var_g.W, N*nu, var);
            }
        }
    } else {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        if(verbose) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
        b1 = HPolySet[0];
        if(!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(ZP, b1, ratios[0], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(ZP, b1, ratios[0], er1, var_g.W, N*nu, var);
        }

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(ZP,HPolySet[i]);
            b2 = HPolySet[i+1];
            if(!window2) {
                var.diameter = 0.6*Rcpp::as<NT>(diam_zono(Rcpp::wrap(ZP.get_mat().transpose()), Rcpp::wrap(A), Rcpp::wrap(zb1.get_vec())));
                std::cout<<"new diameter = "<<var.diameter<<std::endl;
                vol = vol / esti_ratio_interval<RNGType, Point>(zb1, b2, ratios[i], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(zb1, b2, ratios[i], er1, var_g.W, N*nu, var);
            }
        }

        zb1 = ZonoHP(ZP,HPolySet[HPolySet.size()-1]);
        if (!window2) {
            var.diameter = 0.6*Rcpp::as<NT>(diam_zono(Rcpp::wrap(ZP.get_mat().transpose()), Rcpp::wrap(A), Rcpp::wrap(zb1.get_vec())));
            std::cout<<"new diameter = "<<var.diameter<<std::endl;
            vol = vol / esti_ratio_interval<RNGType, Point>(zb1, HP2, ratios[ratios.size() - 1], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(zb1, HP2, ratios[ratios.size() - 1], er1, var_g.W, N*nu, var);
        }
    }

    return vol;

}

#endif
