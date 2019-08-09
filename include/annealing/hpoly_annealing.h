// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALING_H
#define BALL_ANNEALING_H


template <class Polytope, class HPolytope, class VT, typename NT, class Parameters>
void get_hdelta(Polytope &P, HPolytope &HP, VT &Zs_max_gl, NT lb, NT &up_lim, NT &ratio, Parameters &var){

    typedef typename Polytope::PolytopePoint Point;
    typedef typename Polytope::MT MT;
    MT G = P.get_mat().transpose();
    MT A = HP.get_mat();
    int kk = G.cols();
    VT Zs_max = (A*G).cwiseAbs().rowwise().sum();
    Zs_max_gl = Zs_max;
    VT Zs_min = HP.get_vec();
    VT b = HP.get_vec();
    VT b2 = b;
    HPolytope HPiter=HP;

    int n = P.dimension(), m = Zs_max_gl.size();
    int N = 1200;
    Point q(n);
    bool too_few, print = false;
    std::list<Point> randPoints;

    NT l=0.0, u=1.0, med;
    VT  Zmed(m);
    int count =0;
    Parameters variter = var;
    while(true) {

        count++;
        q=Point(n);
        med = (u + l) * 0.5;
        std::cout<<"u = "<<u<<" l = "<<l<<" med = "<<med<<std::endl;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HPiter.set_vec(Zmed);
        variter.che_rad = HPiter.ComputeInnerBall().second;

        randPoints.clear();
        rand_point_generator(HPiter, q, 1200, 10+2*n, randPoints, variter);
        //var.TotSteps = var.TotSteps + NT(1200);

        too_few = false;

        if(check_converg001<Point>(P, randPoints, lb, up_lim, too_few, ratio, 10, 0.2, true, false,var)) {
            HP.set_vec(Zmed);
            return;
        }

        if (too_few) {
            u = med;
        } else {
            l = med;
        }
        if(med>0.9) {
            HP.set_vec(Zmed);
            return;
        }
        if(u-l<0.00000001) {
            std::cout << "fail to find first hpoly... repeat proccess" << std::endl;
            //std::cout<<"origin is in = "<<P.is_in(Point(n))<<std::endl;
            u=1.0;
            l=0.0;
        }
    }
}


template <class Zonotope, class HPolytope, class VT, class PointList, typename NT, class Parameters>
void get_next_zonoball22(Zonotope &Z, std::vector<HPolytope> &HPolySet,
                         HPolytope &HP2, VT Zs_max, VT Zs_min, PointList randPoints,
                        std::vector<NT> &ratios, NT p_value, NT up_lim, int nu, NT alpha, Parameters &var){

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    bool too_few;
    bool print = false;

    NT rad2=0.0;
    NT rad1=0.0, rad;
    NT pnorm, ratio;

    VT Zmed(Zs_max.size());
    NT med, u = 1.0, l = 0.0;

    while (true) {
        med = (u + l) * 0.5;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HP2.set_vec(Zmed);
        too_few = false;

        if(check_converg001<Point>(HP2, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, false,var)){
            HPolySet.push_back(HP2);
            ratios.push_back(ratio);
            return;
        }
        if(too_few) {
            l = med;
        } else {
            u = med;
        }
    }
}

template <class ZonoHP,class Zonotope, class HPolytope, class VT, class Parameters, typename NT>
void get_sequence_of_zonopolys(Zonotope &Z, HPolytope &HP, std::vector<HPolytope> &HPolySet,
                               VT Zs_max, std::vector<NT> &ratios, int Ntot, int nu,
                               NT &p_value, NT up_lim, NT alpha, Parameters &var, Rcpp::Function diam_zono) {

    bool print = var.verbose, too_few=false;
    typedef typename Zonotope::PolytopePoint Point;
    typedef typename Zonotope::MT MT;
    int n = var.n;
    NT ratio, diam_iter;
    std::list<Point> randPoints;
    Point q(n);
    std::cout<<"sampling N points"<<std::endl;
    rand_point_generator(Z, q, Ntot, var.walk_steps, randPoints, var);
    std::cout<<"sampling completed"<<std::endl;
    var.TotSteps = var.TotSteps + NT(Ntot);
    HPolytope HP2 = HP;
    if (check_converg001<Point>(HP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true,var)) {
        ratios.push_back(ratio);
        if(print) std::cout<<"last hpoly and ratio = "<<ratio<<std::endl;
        return;
    }
    if(print) std::cout<<"not the last hpoly"<<std::endl;
    get_next_zonoball22(Z, HPolySet, HP2, Zs_max, HP.get_vec(), randPoints, ratios, p_value, up_lim, nu, alpha, var);
    if(print) std::cout<<"get first hpoly"<<std::endl;

    ZonoHP ZHP2;
    VT Zs_min = HP.get_vec();

    while (true) {

        ZHP2 = ZonoHP(Z,HP2);
        q=Point(n);
        randPoints.clear();

        var.diameter = 0.0;
        for (int j = 0; j < 5; ++j) {
            std::cout<<"computing new diam"<<std::endl;
            diam_iter = Rcpp::as<NT>(
                    diam_zono(Rcpp::wrap(Z.get_mat().transpose()), Rcpp::wrap(HP2.get_mat()*Z.get_mat().transpose()),
                              Rcpp::wrap(HP2.get_vec())));
            // std::cout<<"[annealing] diam_iter = "<<diam_iter<<std::endl;
            if (diam_iter > var.diameter) var.diameter = diam_iter;
        }
        std::cout<<"[annealing] new diameter = "<<var.diameter<<std::endl;

        rand_point_generator(ZHP2, q, Ntot, var.walk_steps, randPoints, var);
        var.TotSteps = var.TotSteps + NT(Ntot);
        if (check_converg001<Point>(HP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true,var)) {
            ratios.push_back(ratio);
            if(print) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
            return;
        }
        get_next_zonoball22(Z, HPolySet, HP2, HP2.get_vec(), Zs_min, randPoints, ratios, p_value, up_lim, nu, alpha,var);
        if(print) std::cout<<"get hpoly"<<std::endl;
    }


}

#endif
