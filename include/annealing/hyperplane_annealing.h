// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#ifndef HYPERPLANE_ANNEALING_H
#define HYPERPLANE_ANNEALING_H





template <class Hyperplane, class Vpolytope, class HPolyBall, class UParameters, typename FT>
void get_next_convex(Vpolytope &P, std::vector<HPolyBall> &PolyBallSet, FT p_value, FT a, UParameters &var, bool &done){

    int n = var.n;
    bool print = var.verbose;
    if (print) std::cout<<"computation of "<<PolyBallSet.size()+1<<" encosing convex body started...\n"<<std::endl;
    HPolyBall Si = PolyBallSet[PolyBallSet.size()-1];
    if (print) std::cout<<"last convex body has = "<<Si.num_of_hyperplanes()<<" hyperplanes"<<std::endl;
    std::list<Point> randPoints;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType &rng2 = var.rng;

    Hyperplane iter_hyp;
    Point q(n);
    if (print) std::cout<<"origin belongs to last convex body: "<<Si.is_in(q)<<std::endl;
    FT rad;
    if (ConvSet.size() ==1) {
        rad = PolyBallSet[0].second().radius();
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere(n, rad));
        }
    } else {
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var);
    }
    q.print();
    if (print) std::cout<<"q after belongs to last convex body: "<<Si.is_in(q)<<std::endl;

    std::list<Point> listIter = randPoints;
    std::list<Point> listIter2;

    int rand_coord;
    //std::vector<int> index(1200, -1);
    int count = 0, count_bef=1200;
    FT totcount = 0, itercount = 0;
    HPolyBall SiIter(n);
    HPolyBall S = Si;
    Point direction(n), cent(n), qbound(n);

    bool added=false;
    bool tested=false, last_hyp = false, last_polyball = false;
    typename  std::list<Point>::iterator pit;
    while (true) {

        //count_bef = count;
        if (print) std::cout<<"points outside new convex = "<<listIter.size()<<std::endl;
        //rand_coord = uidist(rng);
        //index.assign(1200,-1);
        if (print) std::cout<<"q_bef  is in: "<<Si.is_in(q)<<std::endl;
        if (PolyBallSet.size() ==1 && Si.num_of_balls()==1) {
            q = get_point_in_Dsphere(n, rad);
        } else {
            rand_point(Si, q, var);
        }
        if (print) std::cout<<"rand point is in: "<<Si.is_in(q)<<std::endl;
        if (!P.is_in(q)) {
            if (print) std::cout<<"construct ball"<<std::endl;
            //ball_iter = construct_ball(P, q, cent, qbound, direction, var);
            iter_hyp = construct_hyp(P, q, var);
            if (print) std::cout<<"ball constructed"<<std::endl;
            listIter2 = randPoints;
            pit = listIter2.begin();
            while(pit!=listIter2.end()) {
                if (iter_hyp.is_in(*pit)  && Si.is_in(*pit)) {
                    count++;
                    pit=listIter2.erase(pit);
                } else {
                    pit++;
                }
            }
            listIter = listIter2;
            itercount = 0.0;
            //totcount = 0.0;
            std::cout<<"count = "<<count<<std::endl;
            if (count<200) {
                std::cout<<"count = "<<count<<std::endl;
                SiIter = Si;
                SiIter.add_hyperplane(iter_hyp.get_normal_vec(), get_constant());
                if (is_last_hyp(S, SiIter, q, p_value, a, var)) {
                    Si.add_hyperplane(iter_hyp.get_normal_vec(), get_constant());
                    PolyBallSet.push_back(Si);
                    std::cout<<"\nADD CONVEX No. "<<PolyBallSet.size()<<std::endl;
                    return;
                }
                shift_hyperplane(P, randPoints, S, Si, iter_hyp, a, PolyBallSet, last_hyp, last_polyball, var);
                if (last_polyball){
                    done = true;
                    return;
                }
                if (last_hyp) return;
                //ConvSet.push_back(Si);
                //return;
                q = Point(n);
                continue;
            }
            std::cout<<"count_bef - count = "<<count_bef - count<<std::endl;
            if (added && count_bef-count<80) {
                if (is_last_polyball(P, Si, p_value, a, var, false)) {
                    PolyBallSet.push_back(Si);
                    std::cout<<"\nADD LAST CONVEX No. "<<ConvSet.size()<<std::endl;
                    done = true;
                    return;
                } else {
                    added = false;
                }
            }
            q = Point(n);
            Si.add_hyperplane(iter_hyp.get_normal_vec(), get_constant());
            added = true;
            tested = false;
            totcount =0.0;
            count_bef = count;
            count = 0;
        } else {
            itercount += 1.0;
            totcount += 1.0;
            if (totcount>2.0 && !tested) {
                if (is_last_polyball(P, Si, p_value, a, var, false)){
                    PolyBallSet.push_back(Si);
                    done = true;
                    return;
                }
                tested = true;
            }
        }
    }
    //ConvSet.push_back(Si);
}


template <class Hyperplane, class Vpolytope, class HPolyBall, class UParameters, class GParameters, typename FT>
void get_hyperplane_annealing(Vpolytope &VP, std::vector<HPolyBall> &PolyBallSet,
                              FT &t_value, FT &a, UParameters &var1, GParameters &var2) {

    int n = var1.n;
    typedef typename HPolyBall::Polytope Hpolytope;
    std::vector<Ball> vecBalls;
    get_first_conv(VP, vecBalls, var);
    Hpolytope HP(n);
    PolyBallSet.push_back(HPolyBall(HP, vecBalls[0]));

    bool done = false;
    while(!done) {
        get_next_polyball<Hyperplane>(P, PolyBallSet, t_value, a, var1, done);
        if (print) std::cout<<"computation of "<<ConvSet.size()<<" encosing convex body completed.."<<std::endl;
    }


}


#endif
