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

template <class HPolyBall, class VPolytope, class PointList, class Point, class Parameters>
void construct_new_hyp(HPolyBall Si, VPolytope VP, std::vector<HPolyBall> PolyBallSet,
                       PointList &randPoints, Point &q, bool &added_in_set, Parameters var){

    int n = var.n;
    typedef typename Polytope::NT NT;
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    Point center(n), c0(n);
    Point center2(n);
    Point temp(n);
    std::vector <NT> lambdas(P.num_of_vertices());
    //std::cout<<P.is_in(center)<<std::endl;
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    MT Mat = MT::Zero(n,n);
    VT b = VT::Ones(n);
    VT p(n);
    MT V = P.get_mat();
    int count = 0;
    NT z0;
    bool added = false;
    NT min_plus;
    while(!added) {
        min_plus = intersect_line_Vpoly2<NT>(P.get_mat(), center2, v, false, lambdas);


        for (int j = 0; j < P.num_of_vertices(); ++j) {
            if (lambdas[j] > 0.0) {
                Mat.row(count) = V.row(j);
                count++;
            } else {
                p = V.row(j);
            }
        }
        if (count == n) {
            std::cout<<count<<std::endl;
            added = true;
        } else {
            std::cout<<count<<std::endl;
            exit(-1);
        }
    }
    VT a = Mat.colPivHouseholderQr().solve(b);
    if (a.dot(p) > 1.0) {
        a = -a;
        z0 = -1.0;
    }

    HPolyBall SiTemp = Si;
    SiTemp.add_facet(a, z0);

    check_convergence(SiTemp, randPoints, a, done, too_few, var);

    if (!done && !too_few) {
        Si.add_faet(a, z0);
        return;
    }

    if (!too_few && done) {
        Si.add_facet;
        PolyBallSet.push_back(Si);
        added_in_set = true;
        return;
    }

    if (too_few) {
        
    }

}


template <class Hyperplane, class Vpolytope, class HPolyBall, class UParameters, typename FT>
void get_next_convex(Vpolytope &P, std::vector<HPolyBall> &PolyBallSet, FT p_value,
                     FT a, UParameters &var, bool &done){

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
    if (PolyBallSet.size() == 1) {
        rad = PolyBallSet[0].second().radius();
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere(n, rad));
        }
    } else {
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var);
    }
    //q.print();
    //if (print) std::cout<<"q after belongs to last convex body: "<<Si.is_in(q)<<std::endl;

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
    q = Point(n);
    int thre = int(1.0/(1.0-epsilon));

    while (true) {

        do {

            if (count>=thre) {
                check_convergence(Si, VP, epsilon, done, var);
                if (done){
                    return;
                } else {
                    count=-1000;
                }
            }
            count++;
            rand_point(Si, q, var);
        } while (VP.is_in(q)==-1);

        construct_new_hyp(Si, VP, PolyBallSet, randPoints, added_in_set, var);

        if (added_in_set) {
            randPoints.clear();
            q = Point(n);
            rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var);
            count = 0;
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
        get_next_polyball<Hyperplane>(VP, PolyBallSet, t_value, a, var1, done);
        if (print) std::cout<<"computation of "<<ConvSet.size()<<" encosing convex body completed.."<<std::endl;
    }


}


#endif
