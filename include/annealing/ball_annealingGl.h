// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALINGGL_H
#define BALL_ANNEALINGGL_H


template <class Point, class ConvexBody, class PointList, typename NT, class Parameters>
bool check_converg001(ConvexBody &P, PointList &randPoints, NT lb, NT ub, bool &too_few, NT &ratio,
                      int nu, NT alpha, bool precheck, bool lastball, Parameters &var) {

    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int m = randPoints.size()/nu;
    NT T, rs, alpha_check = 0.01, countsIn = 0.0;

    bool print = true;

    int i = 1;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){

        if (P.is_in(*pit)==-1) countsIn += 1.0;
        var.MemLps = var.MemLps + 1.0;
        if (i % m == 0) {
            ratios.push_back(countsIn/m);
            countsIn = 0.0;
            if (ratios.size()>1 && precheck) {
                boost::math::students_t dist(ratios.size() - 1);
                mv = getMeanVariance(ratios);
                ratio = mv.first;
                //std::cout<<"precheck ratio = "<<ratio<<std::endl;
                rs = std::sqrt(mv.second);
                T = rs * (boost::math::quantile(boost::math::complement(dist, alpha_check / 2.0))
                          / std::sqrt(NT(ratios.size())));
                if (ratio + T < lb) {
                    too_few = true;
                    return false;
                } else if (ratio - T > ub) return false;
            }
        }
    }

    //NT alpha = 0.25;
    if(precheck) alpha *= 0.5;
    mv = getMeanVariance(ratios);
    ratio = mv.first;
    //std::cout<<"ratio = "<<ratio<<std::endl;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs*(boost::math::quantile(boost::math::complement(dist, alpha))
              / std::sqrt(NT(nu)));
    if (ratio > lb + T) {
        if (lastball) return true;
        if ((precheck && ratio < ub - T) || (!precheck && ratio < ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;

}


template <class Point, class ball, class PointList, typename NT, class Parameters>
void get_next_zonoball(std::vector<ball> &BallSet, PointList &randPoints, NT rad_min, std::vector<NT> &ratios, NT lb,
                       NT ub, NT alpha, int nu, Parameters &var){

    int n = (*randPoints.begin()).dimension();
    bool too_few;
    NT radmax = 0.0, rad, pnorm, ratio;

    for (typename PointList::iterator rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit) {
        pnorm = (*rpit).squared_length();
        if (pnorm > radmax) radmax = pnorm;
    }
    ball Biter;
    radmax=std::sqrt(radmax);

    while (true) {
        rad = 0.5 * (rad_min + radmax);
        Biter = ball(Point(n), rad * rad);
        too_few = false;

        if (check_converg001<Point>(Biter, randPoints, lb, ub, too_few, ratio, nu, alpha, false, false, var)) {
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return;
        }

        if (too_few) {
            rad_min = rad;
        } else {
            radmax = rad;
        }
    }

}

template <class RNGType,class ball, class Polytope, typename NT, class Parameters>
void get_first_ball(Polytope &P, ball &B0, NT &ratio, NT radius, NT lb, NT ub, NT alpha, NT rmax, Parameters &var){

    typedef typename Polytope::PolytopePoint Point;
    int n = P.dimension();
    bool bisection_int = false, pass = false, too_few = false;
    bool print = true;
    std::list<Point> randPoints;
    Point p(n);

    //std::cout<<"rmax = "<<rmax<<" radius = "<<radius<<std::endl;
    if(rmax>0.0) {
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));
        }
        pass = check_converg001<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false, var);
        if (pass || !too_few) {
            B0 = ball(Point(n), rmax*rmax);
            return;
        }
        bisection_int = true;
    } else {
        rmax = 2 * std::sqrt(NT(n)) * radius;
    }
    NT rad1 = radius;
    //std::cout<<"rmax = "<<rmax<<" rad1 = "<<rad1<<std::endl;

    while(!bisection_int) {
        //std::cout<<"ole!"<<" rmax = "<<rmax<<" rad1 = "<<rad1<<std::endl;

        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));

        if(check_converg001<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false, var)) {
            B0 = ball(Point(n), rmax*rmax);
            return;
        }

        if (too_few) break;
        rad1 = rmax;
        rmax = rmax + 2*std::sqrt(NT(n))*radius;
    }

    NT rad_med;

    while(true) {

        rad_med = 0.5*(rad1+rmax);
        //std::cout<<"rmax = "<<rmax<<" rad1 = "<<rad1<<" rad_med = "<<rad_med<<std::endl;
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rad_med));

        if(check_converg001<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false, var)) {
            B0 = ball(Point(n), rad_med*rad_med);
            //std::cout<<"done, ratio = "<<ratio<<std::endl;
            return;
        }

        if (too_few) {
            rmax = rad_med;
        } else {
            rad1 = rad_med;
        }

    }

}

template <class PolyBall, class RNGType,class ball, class Polytope, class Parameters, typename NT>
void get_sequence_of_polyballs(Polytope &P, std::vector<ball> &BallSet, std::vector<NT> &ratios, int Ntot, int nu,
                               NT lb, NT ub, NT radius, NT alpha, Parameters &var, NT rmax = 0.0, NT ii=0.3) {

    typedef typename Polytope::PolytopePoint Point;
    typedef typename Polytope::MT MT;
    bool print = var.verbose, fail;
    //print = true;
    int n = P.dimension();
    NT ratio, ratio0;
    std::list<Point> randPoints;
    ball B0;
    Point q(n);
    PolyBall zb_it;
    get_first_ball<RNGType>(P, B0, ratio, radius, lb, ub, alpha, rmax, var);
    //std::cout<<"first ball computed"<<std::endl;
    ratio0 = ratio;
    //std::cout<<"is_in = "<<P.is_in(q)<<std::endl;
    rand_point_generator(P, q, Ntot, var.walk_steps, randPoints, var);
    var.TotSteps = var.TotSteps + NT(Ntot);
    //std::cout<<"N ="<<Ntot<<std::endl;

    if (check_converg001<Point>(B0, randPoints, lb, ub, fail, ratio, nu, alpha, false, true, var)) {
        ratios.push_back(ratio);
        BallSet.push_back(B0);
        ratios.push_back(ratio0);
        if(print) std::cout<<"one ball and ratio = "<<ratio<<std::endl;
        return;
    }
    if(print) std::cout<<"not the last ball, ratio = "<<ratio<<std::endl;
    get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, alpha, nu, var);
    if(print) std::cout<<"number of balls = "<<BallSet.size()+1<<std::endl;

    while (true) {
        zb_it = PolyBall(P, BallSet[BallSet.size()-1]);
        var.diameter = ii*2.0*zb_it.radius();
        q=Point(n);
        randPoints.clear();
        //std::cout<<"is_in = "<<zb_it.is_in(q)<<std::endl;
        rand_point_generator(zb_it, q, Ntot, var.walk_steps, randPoints,var);
        var.TotSteps = var.TotSteps + NT(Ntot);
        //std::cout<<"N points sampled from BP"<<std::endl;

        if (check_converg001<Point>(B0, randPoints, lb, ub, fail, ratio, nu, alpha, false, true, var)) {
            ratios.push_back(ratio);
            BallSet.push_back(B0);
            ratios.push_back(ratio0);
            return;
        }
        get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, alpha, nu, var);
        if(print) std::cout<<"number of balls = "<<BallSet.size()+1<<std::endl;
    }
}


#endif
