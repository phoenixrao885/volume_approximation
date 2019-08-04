// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#define VOLESTI_DEBUG
#include "volume.h"
#include "ball_ann_vol.h"
#include "hzono_vol.h"
#include "compute_miniball.h"


template <class Point, class NT, class Polytope>
Rcpp::NumericVector generic_volume(Polytope& P, unsigned int walk_step, double e, Rcpp::Nullable<Rcpp::NumericVector> InnerBall,
                      bool CG, bool BAN, bool hpoly, unsigned int win_len, unsigned int N, double C, double ratio,
                      double frac, double lb, double ub, double p, double alpha, double rmax, unsigned int NN,
                      unsigned int nu, bool win2, bool ball_walk, double delta, bool cdhr, bool rdhr, bool billiard,
                      bool rounding, int type, Rcpp::Nullable<double> diameter, double ii, Rcpp::Function diam_zono)
{

    typedef typename Polytope::VT VT;

    bool rand_only=false,
         birk=false,
         verbose =false;
    unsigned int n_threads=1;

    //unsigned int m;//=A.nrow()-1;
    unsigned int n = P.dimension();//=A.ncol()-1;
    unsigned int rnum = std::pow(e,-2) * 400 * n * std::log(n);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    typedef boost::mt19937    RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::pair<Point,NT> InnerB;
    NT round_val = 1.0, diam, nballs = 0.0;

    if(InnerBall.isNotNull()) { //if it is given as an input
        // store internal point hat is given as input
        Rcpp::NumericVector InnerVec = Rcpp::as<Rcpp::NumericVector>(InnerBall);
        std::vector<NT> temp_p;
        for (unsigned int j=0; j<n; j++){
            temp_p.push_back(InnerVec[j]);
        }
        InnerB.first = Point( n , temp_p.begin() , temp_p.end() );
        // store the radius of the internal ball that is given as input
        InnerB.second = InnerVec[n];
        diam = (diameter.isNotNull()) ? Rcpp::as<NT>(diameter) : 2*InnerB.second;
    }else {
        if (type == 2 && BAN) {
            if (rounding) {
                InnerB.first = P.get_mean_of_vertices();
                InnerB.second = 0.0;
                vars <NT, RNGType> var2(1, n, 1, n_threads, 0.0, e, 0, 0.0, 0, InnerB.second, 2*P.get_max_vert_norm(), rng, urdist, urdist1,
                                        -1, verbose, rand_only, rounding, NN, birk, ball_walk, cdhr, rdhr, billiard,0.0,0.0,0.0);
                std::pair <NT, NT> res_round = rounding_min_ellipsoid(P, InnerB, var2);
                round_val = res_round.first;
                rounding = false;
                InnerB = compute_minball<Point, VT, NT>(P);
                InnerB.second = 0.0;
                //InnerBall.first = Point(n);
                //InnerBall.second = 0.0;
                rmax = P.get_max_vert_norm();
                //std::cout<<"is_in = "<<P.is_in(Point(n))<<std::endl;
                //P.get_mean_of_vertices().print();
            } else {
                InnerB = compute_minball<Point, VT, NT>(P);
                InnerB.second = 0.0;
                rmax = P.get_max_vert_norm();

                //InnerBall.second = 0.0;
                //VP.print();
            }
            //std::cout<<"ii = "<<ii<<std::endl;
            diam = (diameter.isNotNull()) ? Rcpp::as<NT>(diameter) : ii*P.comp_diam();
            //std::cout<<"diameter = "<<diam<<std::endl;
        } else {
            // no internal ball or point is given as input
            InnerB = P.ComputeInnerBall();
            diam = (diameter.isNotNull()) ? Rcpp::as<NT>(diameter) : 2*InnerB.second;
        }

    }

    // initialization
    vars<NT, RNGType> var(rnum,n,walk_step,n_threads,0.0,e,0,0.0,0, InnerB.second, diam, rng,urdist,urdist1,
                          delta,verbose,rand_only,rounding,NN,birk,ball_walk,cdhr,rdhr,billiard, 0.0, 0.0, 0.0);
    NT vol;
    if (CG) {
        vars<NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerB.second, diam, rng,
                               urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk,
                               cdhr,rdhr,billiard, 0.0, 0.0, 0.0);
        vars_g<NT, RNGType> var1(n, walk_step, N, win_len, 1, e, InnerB.second, rng, C, frac, ratio, delta, false, verbose,
                                 rand_only, rounding, false, birk, ball_walk, cdhr, rdhr, 0.0, 0.0);
        vol = volume_gaussian_annealing(P, var1, var2, InnerB, nballs);
        var.BoundCalls =  var1.BoundCalls;
        var.TotSteps = var1.TotSteps;
    } else if (BAN) {
        vars_ban <NT> var_ban(lb, ub, p, rmax, alpha, win_len, NN, nu, win2);
        if (!hpoly) {
            //vars_ban <NT> var_ban(0.1, 0.15, 0.75, 0.0, 0.0, 0.0, 2 * n * n + 250, 220 + (n * n) / 10, 10, false);
            //std::cout<<"ii = "<<ii<<std::endl;
            vol = volesti_ball_ann(P, var, var_ban, InnerB, nballs, ii);
        } else {
            vars_g<NT, RNGType> varg(n, 1, N, 6*n*n+500, 1, e, InnerB.second, rng, C, frac, ratio, delta, false, verbose,
                                     rand_only, false, false, birk, false, true, false, 0.0, 0.0);
            vol = vol_hzono<HPolytope<Point> > (P, var, var_ban, varg, InnerB,nballs,diam_zono);
        }
    }else {
        vol = volume(P, var, InnerB);
    }

    Rcpp::NumericVector res(5);
    res[0] = vol*round_val;
    res[1] = nballs;
    res[2] = var.BoundCalls;
    res[3] = var.MemLps;
    res[4] = var.TotSteps;

    return res;
    //return round_val * vol;
}

//' The main function for volume approximation of a convex Polytope (H-polytope, V-polytope or a zonotope)
//'
//' For the volume approximation can be used two algorithms. Either SequenceOfBalls or CoolingGaussian. A H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}. A V-polytope is defined as the convex hull of \eqn{m} \eqn{d}-dimensional points which correspond to the vertices of P. A zonotope is desrcibed by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param walk_step Optional. The number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for SequenceOfBalls and \eqn{1} for CoolingGaussian.
//' @param error Optional. Declare the upper bound for the approximation error. The default value is \eqn{1} for SequenceOfBalls and \eqn{0.1} for CoolingGaussian.
//' @param InnerBall Optional. A \eqn{d+1} vector that contains an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an inscribed ball.
//' @param Algo Optional. A string that declares which algorithm to use: a) \code{'SoB'} for SequenceOfBalls or b) \code{'CG'} for CoolingGaussian or c) \code{'BAN'} for CoolingBodies.
//' @param WalkType Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run or c) \code{'BW'} for Ball Walk. The default walk is \code{'CDHR'}.
//' @param rounding Optional. A boolean parameter for rounding. The default value is \code{FALSE}.
//' @param Parameters Optional. A list for the parameters of the algorithms:
//' \itemize{
//' \item{\code{Window} }{ The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
//'  \item{\code{C} }{ A constant for the lower bound of \eqn{variance/mean^2} in schedule annealing of CG algorithm. The default value is \eqn{2}.}
//'  \item{\code{M} }{ The number of points we sample in each step of schedule annealing in CG algorithm. The default value is \eqn{500C + dimension^2 / 2}.}
//'  \item{\code{ratio} }{ Parameter of schedule annealing of CG algorithm, larger ratio means larger steps in schedule annealing. The default value is \eqn{1 - 1/dimension}.}
//'  \item{\code{frac} }{ The fraction of the total error to spend in the first gaussian in CG algorithm. The default value is \eqn{0.1}.}
//'  \item{\code{BW_rad} }{ The radius for the ball walk. The default value is \eqn{4r/dimension}, where \eqn{r} is the radius of the inscribed ball of the polytope.}
//'  \item{\code{ub} }{ The lower bound for the ratios in MMC in BAN algorithm. The default value is \eqn{0.1}.}
//'  \item{\code{lb} }{ The upper bound for the ratios in MMC in BAN algorithm. The default value is \eqn{0.15}.}
//'  \item{\code{N} }{ An integer that controls the number of points \eqn{\nu N} generated in each convex body in annealing schedule.}
//'  \item{\code{nu} }{ The degrees of freedom for the t-student distribution in t-tests in BAN algorithm. The default value is \eqn{10}.}
//'  \item{\code{alpha} }{ The significance level for the t-tests in BAN algorithm. The default values is 0.2.}
//'  \item{\code{prob} }{ The probability is used for the empirical confidence interval in ratio estimation of BAN algorithm. The default value is \eqn{0.75}.}
//'  \item{\code{hpoly} }{ A boolean parameter to use H-polytopes in MMC of BAN algorithm. The default value is \code{FALSE}.}
//'  \item{\code{minmaxW} }{ A boolean parameter to use the sliding window with a minmax values stopping criterion.}
//'
//' }
//'
//' @references \cite{I.Z.Emiris and V. Fisikopoulos,
//' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2014.},
//' @references \cite{B. Cousins and S. Vempala, \dQuote{A practical volume algorithm,} \emph{Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society,} 2015.}
//'
//'
//' @return The approximation of the volume of a convex polytope.
//' @examples
//' # calling SOB algorithm for a H-polytope (2d unit simplex)
//' P = GenSimplex(2,'H')
//' vol = volume(P)
//'
//' # calling CG algorithm for a V-polytope (3d simplex)
//' P = GenSimplex(2,'V')
//' vol = volume(P, Algo = "CG")
//'
//' # calling CG algorithm for a 2-dimensional zonotope defined as the Minkowski sum of 4 segments
//' Z = GenZonotope(2, 4)
//' vol = volume(Z, WalkType = "RDHR", walk_step = 5)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector volume (Rcpp::Reference P, Rcpp::Function diam_zono, Rcpp::Nullable<unsigned int> walk_step = R_NilValue,
                Rcpp::Nullable<double> error = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> InnerBall = R_NilValue,
                Rcpp::Nullable<std::string> Algo = R_NilValue,
                Rcpp::Nullable<std::string> WalkType = R_NilValue,
                Rcpp::Nullable<bool> rounding = R_NilValue,
                Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue,
                Rcpp::Nullable<double> diameter = R_NilValue,
                Rcpp::Nullable<double> iii = R_NilValue) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef IntersectionOfVpoly <Vpolytope> InterVP;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    unsigned int n = P.field("dimension"), walkL;

    bool CG, BAN, cdhr = true, rdhr = false, ball_walk = false, billiard = false, round, win2 = false, hpoly = false;
    unsigned int win_len = 4 * n * n + 500, N = 500 * 2 + n * n / 2, NN = 120 + (n*n)/10, nu = 10;

    double C = 2.0, ratio = 1.0 - 1.0 / (NT(n)), frac = 0.1, e, delta = -1.0,
            lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2, ii=0.4;

    if (!rounding.isNotNull()) {
        round = false;
    } else {
        round = Rcpp::as<bool>(rounding);
    }

    if (!WalkType.isNotNull() || Rcpp::as<std::string>(WalkType).compare(std::string("CDHR")) == 0) {
        cdhr = true;
        rdhr = false;
        ball_walk = false;
    } else if (Rcpp::as<std::string>(WalkType).compare(std::string("RDHR")) == 0) {
        cdhr = false;
        rdhr = true;
        ball_walk = false;
    } else if (Rcpp::as<std::string>(WalkType).compare(std::string("BilW"))==0) {
        cdhr = false;
        rdhr = false;
        ball_walk = false;
        billiard = true;
    } else if (Rcpp::as<std::string>(WalkType).compare(std::string("BW"))==0) {
        cdhr = false;
        rdhr = false;
        ball_walk = true;
    } else {
        throw Rcpp::exception("Unknown walk type!");
    }

    if (!Algo.isNotNull() || Rcpp::as<std::string>(Algo).compare(std::string("SOB")) == 0) {

        CG = false;
        BAN = false;

        if (!walk_step.isNotNull()) {
            walkL = 10 + n / 10;
        } else {
            walkL = Rcpp::as<unsigned int>(walk_step);
        }

        if (!error.isNotNull()) {
            e = 1.0;
        } else {
            e = Rcpp::as<NT>(error);
        }

    } else if (Rcpp::as<std::string>(Algo).compare(std::string("CG")) == 0) {

        CG = true;
        BAN = false;

        if (!error.isNotNull()) {
            e = 0.1;
        } else {
            e = Rcpp::as<NT>(error);
        }

        if (!walk_step.isNotNull()) {
            walkL = 1;
        } else {
            walkL = Rcpp::as<int>(walk_step);
        }

    } else if (Rcpp::as<std::string>(Algo).compare(std::string("BAN")) == 0) {
        CG = false;
        BAN = true;
        if (!error.isNotNull()) {
            e = 0.1;
        } else {
            e = Rcpp::as<NT>(error);
        }
        if (!walk_step.isNotNull()) {
            walkL = 1;
        } else {
            walkL = Rcpp::as<int>(walk_step);
        }

    }else {
        throw Rcpp::exception("Unknown method!");
    }
    if (iii.isNotNull()) {
        ii = Rcpp::as<NT>(iii);
    }

    if(Parameters.isNotNull()) {

        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("BW_rad")) {
            delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["BW_rad"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("C")) {
            C = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["C"]);
            N = 500 * ((int) C) + n * n / 2;
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("M")) {
            N = Rcpp::as<int>(Rcpp::as<Rcpp::List>(Parameters)["M"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("Window")) {
            win_len = Rcpp::as<int>(Rcpp::as<Rcpp::List>(Parameters)["Window"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("frac")) {
            frac = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["frac"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("ratio")) {
            ratio = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["ratio"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("hpoly")) {
            hpoly = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(Parameters)["hpoly"]);
            if ((hpoly && !BAN) || (Rcpp::as<int>(P.field("type"))!=3 && BAN && hpoly)) Rf_warning("flag 'hpoly' can be used to only"
                                                                                    " in MMC of BAN algorithm for"
                                                                                    " zonotopes.");
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("lb")) {
            lb = Rcpp::as<double>(Rcpp::as<Rcpp::List>(Parameters)["lb"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("ub")) {
            ub = Rcpp::as<double>(Rcpp::as<Rcpp::List>(Parameters)["ub"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("nu")) {
            nu = Rcpp::as<int>(Rcpp::as<Rcpp::List>(Parameters)["nu"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("N")) {
            NN = Rcpp::as<int>(Rcpp::as<Rcpp::List>(Parameters)["N"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("minmaxW")) {
            win2 = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(Parameters)["minmaxW"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("prob")) {
            p = Rcpp::as<double>(Rcpp::as<Rcpp::List>(Parameters)["prob"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("alpha")) {
            alpha = Rcpp::as<double>(Rcpp::as<Rcpp::List>(Parameters)["alpha"]);
        }
    }

    int type = P.field("type");
    switch(type) {
        case 1: {
            // Hpolytope
            Hpolytope HP;
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));

            return generic_volume<Point, NT>(HP, walkL, e, InnerBall, CG, BAN, hpoly, win_len, N, C, ratio, frac, lb,
                                             ub, p, alpha, rmax, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard,
                                             round, type, diameter, ii, diam_zono);
        }
        case 2: {
            // Vpolytope
            Vpolytope VP;
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));

            return generic_volume<Point, NT>(VP, walkL, e, InnerBall, CG, BAN, hpoly, win_len, N, C, ratio, frac, lb,
                                             ub, p, alpha, rmax, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard,
                                             round, type, diameter, ii, diam_zono);
        }
        case 3: {
            // Zonotope
            zonotope ZP;
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));

            return generic_volume<Point, NT>(ZP, walkL, e, InnerBall, CG, BAN, hpoly, win_len, N, C, ratio, frac, lb,
                                             ub, p, alpha, rmax, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard,
                                             round, type, diameter, ii, diam_zono);
        }
        case 4: {
            // Intersection of two V-polytopes
            Vpolytope VP1;
            Vpolytope VP2;
            InterVP VPcVP;
            VP1.init(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
            VP2.init(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
            VPcVP.init(VP1, VP2);
            Rcpp::NumericVector InnerVec(n+1);
            bool empty;
            std::pair<Point, NT> InnerB = VPcVP.getInnerPoint_rad(empty);
            if (empty) {
                Rf_warning("Empty set");
                return 0;
            }
            if(!InnerBall.isNotNull()) {
                for (int i = 0; i < n; ++i) {
                    InnerVec[i] = InnerB.first[i];
                }
                InnerVec[n] = InnerB.second;
            } else {

                InnerVec = Rcpp::as<Rcpp::NumericVector>(InnerBall);

            }

            return generic_volume<Point, NT>(VPcVP, walkL, e, InnerVec, CG, BAN, hpoly, win_len, N, C, ratio, frac, lb,
                                             ub, p, alpha, rmax, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard,
                                             round, type, diameter, ii, diam_zono);
        }
    }

    return Rcpp::NumericVector(4);
}
