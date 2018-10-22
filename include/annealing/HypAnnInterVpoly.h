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

#ifndef HYPANNINTERVPOLY_H
#define HYPANNINTERVPOLY_H

template <class Polytope, class Ball, class Point, class Parameters>
void enclosing_ball2(Polytope &P, Ball &B0, Point &xc, Parameters &var){

    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Polytope::NT NT;
    //typedef typename Polytope::PolytopePoint Point;
    unsigned int n = P.dimension();
    std::vector<NT> vec(n,0.0);
    NT rad=0.0, nr;
    VT c_e(n);
    for (unsigned int k = 0; k < n; ++k) {
        c_e(k) = xc[k];
    }
    P.shift(c_e);
    MT V1 = P.get_mat1();
    MT V2 = P.get_mat2();
    unsigned int k1 = V1.rows();
    unsigned int k2 = V2.rows();
    Point temp(n);

    for (int i = 0; i < k1; ++i) {
        for (int j = 0; j < n; ++j) {
            vec[j] = V1(i,j);
        }
        temp = Point(n, vec.begin(), vec.end());
        nr = std::sqrt(temp.squared_length());
        if ( nr > rad) {
            rad = nr;
        }
    }

    for (int i = 0; i < k2; ++i) {
        for (int j = 0; j < n; ++j) {
            vec[j] = V2(i,j);
        }
        temp = Point(n, vec.begin(), vec.end());
        nr = std::sqrt(temp.squared_length());
        if ( nr > rad) {
            rad = nr;
        }
    }

    B0 = Ball(Point(n), rad*rad);

}


template <class Point, class HPolyBall, class PointList, typename NT, class Parameters>
void check_converg3(HPolyBall &Si, PointList &randPoints, NT p_test, bool &done, bool &too_few, Parameters &var, bool print) {

    //typedef typename HPolyBall::Point Point;
    std::vector<NT> ratios;
    NT countsIn = 0.0;

    int i = 1;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (Si.is_in(*pit)) {
            countsIn += 1.0;
        }
        if (i % 120 == 0) {
            if (print) std::cout<<"ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countsIn/120.0);
            countsIn = 0.0;
        }
    }

    std::pair<NT,NT> mv = getMeanVariance(ratios);
    NT t_value = 0.700;
    NT p_mval = mv.first;
    NT p_varval = mv.second;
    int ni = ratios.size();
    //NT p_test = a;

    if (print) std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
    if (p_mval > p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))) {
        if (p_mval < (p_test + 0.05) + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))) {
            done= true;
        }
    } else {
        too_few = true;
    }

}


template <class HPolyBall, class VPolytope, class PointList, class VT, class Point, typename NT, class Parameters>
void reconstruct_hyp2(HPolyBall &Si, VPolytope &VP, std::vector<HPolyBall> &PolyBallSet,
                     PointList &randPoints, VT &c, NT &z0, Point &v, NT &p_value, Parameters &var){

    //typedef typename VPolytope::PolytopePoint Point;
    unsigned int n = var.n;
    std::pair<NT,NT> bpair = PolyBallSet[0].line_intersect(Point(n),v);
    Point q = bpair.first*v;
    NT zmax = 0.0, zmed;
    bool doneIter, too_few;
    Parameters var2 = var;
    var2.coordinate = true;

    typename std::vector<NT>::iterator pit = q.iter_begin();
    int i = 0;
    for ( ;  pit!=q.iter_end(); ++pit, ++i) {
        zmax += c(i) * (*pit);
    }
    zmax= 10*zmax;
    std::cout<<"z0 = "<<z0<<" zmax = "<<zmax<<std::endl;

    /*NT sum;
    int countin = 0, countin2 = 0;
    for(typename std::list<Point>::iterator prit=randPoints.begin(); prit!=randPoints.end(); ++prit) {
        pit = (*prit).iter_begin();
        int i = 0;
        sum = 0.0;
        for (; pit != (*prit).iter_end(); ++pit, ++i) {
            sum += c(i) * (*pit);
        }
        if (sum<zmax) {
            countin++;
        }
        if (sum<z0) {
            countin2++;
        }
    }*/
    //std::cout<<"inside max halfspace = "<<countin<<std::endl;
    //std::cout<<"inside min halfspace = "<<countin2<<" p_value = "<<p_value<<std::endl;

    NT z00 = z0;
    HPolyBall SiIter;
    while (true) {

        while (true) {

            zmed = (z0 + zmax) / 2.0;
            SiIter = Si;
            doneIter = false;
            too_few = false;

            SiIter.add_facet(c, zmed);
            check_converg3<Point>(SiIter, randPoints, p_value, doneIter, too_few, var, false);

            //std::cout<<"2nd loop.... "<<"doneIter = "<<doneIter<< " too_few = "<<too_few<<std::endl;
            if (doneIter) {
                PolyBallSet.push_back(SiIter);
                randPoints.clear();
                q = Point(n);
                rand_point_generator(SiIter, q, 1200, var.walk_steps, randPoints, var2);
                break;
            }

            if (!doneIter && !too_few) {
                zmax = zmed;
            }

            if (too_few) {
                z0 = zmed;
            }

        }

        SiIter = Si;
        SiIter.add_facet(c, z00);
        doneIter = false;
        too_few = false;
        //std::cout<<randPoints.size()<<std::endl;
        check_converg3<Point>(SiIter, randPoints, p_value, doneIter, too_few, var, false);
        //std::cout<<"doneIter = "<<doneIter<< " too_few = "<<too_few<<std::endl;

        if(doneIter) {
            PolyBallSet.push_back(SiIter);
            Si = SiIter;
            randPoints.clear();
            q = Point(n);
            rand_point_generator(SiIter, q, 1200, var.walk_steps, randPoints, var2);
            return;
        }

        if (!too_few) {
            Si.add_facet(c, z00);
            //randPoints.clear();
            //q = Point(n);
            //rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var2);
            return;
        }

        if (too_few) {
            zmax = zmed;
            z0 = z00;
        }

    }

}


template <class HPolyBall, class VPolytope, class PointList, class Point, typename NT, class Parameters>
void construct_new_hyp2(HPolyBall &Si, VPolytope &VP, std::vector<HPolyBall> &PolyBallSet,
                       PointList &randPoints, Point &q, NT p_value, Parameters &var){

    unsigned int n = var.n;
    //typedef typename Polytope::NT NT;
    typedef typename VPolytope::VT VT;
    typedef typename VPolytope::MT MT;
    Parameters var2 = var;
    var2.coordinate = true;
    Point center(n), c0(n);
    Point center2(n);
    Point temp(n);
    std::vector <NT> lambdas1(VP.first().num_of_vertices());
    std::vector <NT> lambdas2(VP.second().num_of_vertices());
    //std::cout<<P.is_in(center)<<std::endl;
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    MT Mat = MT::Zero(n,n);
    VT b = VT::Ones(n);
    VT p(n);
    MT V1 = VP.get_mat1();
    MT V2 = VP.get_mat2();
    int count = 0;
    NT z0 = 1.0;
    bool added = false, done = false, too_few = false;
    NT min_plus1, min_plus2;

    min_plus1 = intersect_line_Vpoly2<NT>(VP.get_mat1(), center2, v, false, lambdas1);
    min_plus2 = intersect_line_Vpoly2<NT>(VP.get_mat2(), center2, v, false, lambdas2);

    if (min_plus1 < min_plus2) {
        for (int j = 0; j < VP.first().num_of_vertices(); ++j) {
            if (lambdas1[j] > 0.0) {
                Mat.row(count) = V1.row(j);
                count++;
            } else {
                p = V1.row(j);
            }
        }
        if (count != n) {
            std::cout<<count<<std::endl;
            exit(-1);
        }
    } else {
        for (int j = 0; j < VP.second().num_of_vertices(); ++j) {
            if (lambdas2[j] > 0.0) {
                Mat.row(count) = V2.row(j);
                count++;
            } else {
                p = V2.row(j);
            }
        }
        if (count != n) {
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

    check_converg3<Point>(SiTemp, randPoints, p_value, done, too_few, var, false);

    if (!done && !too_few) {
        Si.add_facet(a, z0);
        //randPoints.clear();
        //q = Point(n);
        //rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var2);
        return;
    }

    if (done) {
        Si.add_facet(a, z0);
        PolyBallSet.push_back(Si);
        //added_in_set = true;
        randPoints.clear();
        q = Point(n);
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var2);
        return;
    }

    if (too_few) {
        std::cout<<"reconstructing hyperplanes.."<<std::endl;
        reconstruct_hyp2(Si, VP, PolyBallSet, randPoints, a, z0, v, p_value, var);
        //added_in_set = true;
    }

}


template <class Vpolytope, class HPolyBall, class Parameters, typename NT>
void get_seq_convex2(Vpolytope &VP, std::vector<HPolyBall> &PolyBallSet, NT p_value, NT epsilon,
                    Parameters &var, bool &done){

    typedef typename Vpolytope::PolytopePoint Point;
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n;
    Parameters var2 = var;
    var2.coordinate = true;
    bool print = var.verbose;
    //if (print) std::cout<<"computation of "<<PolyBallSet.size()+1<<" encosing convex body started...\n"<<std::endl;
    HPolyBall Si = PolyBallSet[PolyBallSet.size()-1];
    //if (print) std::cout<<"last convex body has = "<<Si.num_of_hyperplanes()<<" hyperplanes"<<std::endl;
    std::list<Point> randPoints;
    std::list<Point> SirandPoints;
    Point q(n);
    //if (print) std::cout<<"origin belongs to last convex body: "<<Si.is_in(q)<<std::endl;
    NT rad;
    if (PolyBallSet.size() == 1) {
        rad = PolyBallSet[0].second().radius();
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere<RNGType, Point> (n, rad));
        }
    } else {
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var);
    }
    //q.print();
    //if (print) std::cout<<"q after belongs to last convex body: "<<Si.is_in(q)<<std::endl;

    //std::list<Point> listIter = randPoints;
    //std::list<Point> listIter2;

    //int rand_coord;
    //std::vector<int> index(1200, -1);
    int count = 0;//, count_bef=1200;
    //NT totcount = 0, itercount = 0;
    //HPolyBall SiIter(n);
    //HPolyBall S = Si;
    // Point direction(n), cent(n), qbound(n);

    //bool added = false;
    bool too_few = false;
    //typename  std::list<Point>::iterator pit;
    int thre = int(3.0/(1.0-epsilon));

    while (true) {

        q = Point(n);
        do {

            //std::cout<<"count = "<<count<<" epsilon = "<<epsilon<<std::endl;
            if (count>=thre) {
                if (Si.num_of_hyperplanes() > PolyBallSet[PolyBallSet.size()-1].num_of_hyperplanes()) {
                    q = Point(n);
                    SirandPoints.clear();
                    rand_point_generator(Si, q, 1200, var.walk_steps, SirandPoints, var2);
                    check_converg3<Point>(VP, SirandPoints, epsilon, done, too_few, var, false);
                } else {
                    check_converg3<Point>(VP, randPoints, epsilon, done, too_few, var, false);
                }
                //std::cout<<"VP convergence.... "<<"doneIter = "<<done<< " too_few = "<<too_few<<std::endl;
                if (done){
                    if (Si.num_of_hyperplanes() > PolyBallSet[PolyBallSet.size()-1].num_of_hyperplanes()) {
                        PolyBallSet.push_back(Si);
                    }
                    std::cout<<"number of hyps in Si = "<<Si.num_of_hyperplanes()<<std::endl;
                    return;

                } else if(!too_few) {
                    if (Si.num_of_hyperplanes() > PolyBallSet[PolyBallSet.size()-1].num_of_hyperplanes()) {
                        PolyBallSet.push_back(Si);
                    }
                    std::cout<<"number of hyps in Si = "<<Si.num_of_hyperplanes()<<std::endl;
                    done = true;
                    return;
                } else {
                    count=-1000;
                    too_few = false;
                }
            }
            count++;
            rand_point(Si, q, var);
        } while (VP.is_in(q)==-1);

        //std::cout<<"constructing new hyperplane.."<<std::endl;
        construct_new_hyp2(Si, VP, PolyBallSet, randPoints, q, p_value, var);
        std::cout<<"number of hyps = "<<Si.num_of_hyperplanes()<<std::endl;
        count = 0;
        too_few = false;

    }
    //ConvSet.push_back(Si);
}


template <class Vpolytope, class HPolyBall, class Parameters, class Point, typename NT>
void get_hyperplane_annealing2(Vpolytope &VP, std::vector<HPolyBall> &PolyBallSet,
                              NT &p_value, NT &a, NT epsilon, Point &xc, Parameters &var) {

    unsigned int n = var.n;
    typedef typename HPolyBall::HPolytope Hpolytope;
    typedef typename HPolyBall::ball Ball;
    //typedef typename Vpolytope::PolytopePoint Point;
    std::vector<Ball> vecBalls;
    //get_first_conv(VP, vecBalls, var);
    Ball B0;
    Point center(n);
    enclosing_ball2(VP, B0, xc, var);
    Hpolytope HP(n);
    PolyBallSet.push_back(HPolyBall(HP, B0));

    bool done = false;
    get_seq_convex2(VP, PolyBallSet, p_value, epsilon, var, done);

}


#endif
