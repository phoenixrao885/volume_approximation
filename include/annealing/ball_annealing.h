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

#ifndef BALL_ANNEALING_H
#define BALL_ANNEALING_H


template <typename FT>
bool check_t_test(std::vector<FT> ratios, FT a, FT p_test){

    bool passed = false;
    std::pair<FT,FT> mv = getMeanVariance(ratios);
    FT t_value = 0.700;
    FT p_mval = mv.first;
    FT p_varval = mv.second;
    int ni = ratios.size();
    //std::cout<<"mean = "<<mv.first<<" var = "<<mv.second<<std::endl;
    if(mv.first >=190.0/1200.0 && false) {
        return true;
    }

    //std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(FT(ni)))<<std::endl;
    if (p_mval > p_test + t_value*(ni-1)*(p_varval/std::sqrt(FT(ni)))) {
        passed = true;
    }
    return passed;
    //return false;
}


template <class Point, class Interballs, class PointList, typename NT, class Parameters>
void check_converg2(Interballs &Si, PointList &randPoints, NT p_test, bool &done, bool &too_few, Parameters &var, bool print) {

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

//construct_ball(Si, VP, ConvSet, randPoints, q, p_value, var);
template <class Interballs, class VPolytope, class Point, class PointList, typename NT, class Parameters>
void construct_ball(Interballs &Si, VPolytope &VP, std::vector<Interballs> &ConvSet,
                    PointList &randPoints, Point &q, NT &p_value, Parameters &var) {

    int n = var.n;
    //typedef typename VPolytope::NT NT;
    typedef typename VPolytope::VT VT;
    typedef typename VPolytope::MT MT;
    typedef typename Interballs::cball Ball;
    Parameters var2 = var;
    var2.coordinate = true;
    Point center(n), c0(n);
    Point center2(n);
    Point temp(n);
    std::vector <NT> lambdas(VP.num_of_vertices());
    //std::cout<<P.is_in(center)<<std::endl;
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    MT Mat = MT::Zero(n,n);
    VT b = VT::Ones(n);
    VT p(n);
    MT V = VP.get_mat();
    int count = 0;
    NT z0;
    bool added = false;
    NT min_plus;
    while(!added) {
        min_plus = intersect_line_Vpoly2<NT>(VP.get_mat(), center2, v, false, lambdas);


        for (int j = 0; j < VP.num_of_vertices(); ++j) {
            if (lambdas[j] > 0.0) {
                Mat.row(count) = V.row(j);
                count++;
            } else {
                p = V.row(j);
            }
        }
        if (count == n) {
            //std::cout<<count<<std::endl;
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

    //q = v * (min_plus + 0.0001);
    q = v * (min_plus + 0.001);
    typename std::vector<NT>::iterator tempit = temp.iter_begin();
    int i = 0;
    for ( ; tempit!=temp.iter_end(); ++tempit, ++i) {
        *tempit = a(i);
    }
    temp = temp * (1.0 / std::sqrt(temp.squared_length()));
    temp = temp * (min_plus * 50.0);
    int counter = 0;
    NT const tol = 0.000001;
    //temp = q;
    while (true) {
        counter++;
        c0 = center;

        center = center - temp;
        //std::cout<<std::sqrt( (center - q).squared_length())<<std::endl;
        if (VP.is_in_ball(center, std::sqrt( (center - q).squared_length()) )) {
            break;
        }
    }
  //  std::cout<<"first ball found, counter = "<<counter<<std::endl;
    Point midpoint(n);
    counter = 0;
   // std::cout<<"bisection method to find closest center\n";
    while ( std::sqrt((center - c0).squared_length())>tol ) {
        counter++;
        midpoint = (center + c0) * 0.5;
        if (VP.is_in_ball(midpoint, std::sqrt( (midpoint - q).squared_length()) )) {
            center = midpoint;
        } else {
            c0 = midpoint;
        }
    }
  //  std::cout<<"closest center found, counter = "<<counter<<std::endl;
   // center.print();
   // std::cout<<"radius = "<<std::sqrt((center - q).squared_length())<<std::endl;

    Ball B0(center,  (center - q).squared_length() );
    Interballs SiTemp = Si;
    SiTemp.add_ball(B0);
    NT rerad;
    bool done = false, too_few = false;

    check_converg2<Point>(SiTemp, randPoints, p_value, done, too_few, var, false);
   // std::cout<<"initial ball.... "<<"done = "<<done<< " too_few = "<<too_few<<std::endl;

    if (!done && !too_few) {
        Si.add_ball(B0);
        //randPoints.clear();
        //q = Point(n);
        //rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var2);
        return;
    }

    if (done) {
        Si.add_ball(B0);
        ConvSet.push_back(Si);
        //added_in_set = true;
        randPoints.clear();
        q = Point(n);
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var2);
        return;
    }

    if (too_few) {
       // std::cout<<"reconstructing balls.."<<std::endl;
        rerad = B0.radius();
        reconstruct_ball(Si, VP, ConvSet, randPoints, center, rerad, p_value, var);
        //added_in_set = true;
    }

    //return Ball(center,  (center - q).squared_length() );
}


//reconstruct_ball(Si, VP, ConvSet, randPoints, center, q, v, p_value, var);
template <class VPolytope, class PointList, class Interballs, class Point, typename NT, class Parameters>
void reconstruct_ball(Interballs Si, VPolytope &VP, std::vector<Interballs> &Convset,
                      PointList &randPoints, Point &center, NT &radius, NT p_value,
                      Parameters var){

    int n = var.n;
    typedef typename Interballs::cball Ball;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;
    Parameters var2 = var;
    var2.coordinate = true;

    typename  std::list<Point>::iterator rpit = randPoints.begin();
    int i =0;
    NT maxrad = 0.0, temprad;
    for ( ; rpit!=randPoints.end(); ++rpit, ++i) {
        temprad = std::sqrt((center- (*rpit)).squared_length());
        if (temprad>maxrad) maxrad = temprad;
    }

    NT rad1 = radius, rad2 = maxrad, midrad;
    //std::cout<<"rad1 = "<<rad1<<" rad2 = "<<rad2<<std::endl;
    Interballs SiIter;
    Ball ball_iter;
    bool doneIter, too_few;
    Point q(n);

    while (true) {
        //firstit = true;
        while (true) {

            midrad = (rad1 + rad2) / 2.0;
            ball_iter = Ball(center, midrad * midrad);
            doneIter = false;
            too_few = false;

            SiIter.add_ball(ball_iter);
            check_converg2<Point>(SiIter, randPoints, p_value, doneIter, too_few, var, false);

            //std::cout<<"2nd loop.... "<<"doneIter = "<<doneIter<< " too_few = "<<too_few<<std::endl;
            if (doneIter) {
                Convset.push_back(SiIter);
                randPoints.clear();
                q = Point(n);
                rand_point_generator(SiIter, q, 1200, var.walk_steps, randPoints, var2);
                break;
            }

            if (!doneIter && !too_few) {
                rad2 = midrad;
            }

            if (too_few) {
                rad1 = midrad;
            }
        }

        SiIter = Si;
        ball_iter = Ball(center, radius*radius);
        SiIter.add_ball(ball_iter);
        doneIter = false;
        too_few = false;
        //std::cout<<randPoints.size()<<std::endl;
        check_converg2<Point>(SiIter, randPoints, p_value, doneIter, too_few, var, false);
        //std::cout<<"doneIter = "<<doneIter<< " too_few = "<<too_few<<std::endl;

        if (doneIter) {
            Convset.push_back(SiIter);
            Si = SiIter;
            randPoints.clear();
            q = Point(n);
            rand_point_generator(SiIter, q, 1200, var.walk_steps, randPoints, var2);
            return;
        }

        if (!too_few) {
            Si.add_ball(ball_iter);
            //randPoints.clear();
            //q = Point(n);
            //rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var2);
            return;
        }

        if (too_few) {
            rad2 = midrad;
            rad1 = radius;
        }
    }
}



template <class Polytope, class Interballs, class Ball, class Parameters>
void get_first_conv(Polytope &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls, Parameters &var) {

    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Polytope::NT NT;
    typedef typename Polytope::PolytopePoint Point;
    unsigned int n = P.dimension();
    std::cout<<"dimension = "<<n<<std::endl;
    std::vector<NT> vec(n,0.0);

    Point xc(n);

    MT V = P.get_mat();
    P.print();
    int k = V.rows();
    Point temp;

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            vec[j] = V(i,j);
        }
        temp = Point(n, vec.begin(), vec.end());
        xc = xc + temp;
    }
    xc = xc * (1.0/NT(k));

    NT rad = 0.0;
    NT nr;
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            vec[j] = V(i,j);
        }
        temp = Point(n, vec.begin(), vec.end());
        nr = std::sqrt((xc - temp).squared_length());
        if ( nr > rad) {
            rad = nr;
        }
    }
    rad = rad*1.000001;
    VT c_e(n);
    for(int i=0; i<n; i++){
        c_e(i)=xc[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);
    //std::cout<<"radius of minim ball = "<<rad<<std::endl;
    xc.print();
    //Point xc = Point(n);
    //std::vector<Ball> S0;
    vecBalls.push_back(Ball(Point(n), rad*rad));
    //ConvSet.push_back(Interballs(P.dimension(), vecBalls));

}
/*
template <class T>
void get_first_conv(T &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls) {
    std::pair<Point, NT> res = compute_minball(P);
    int n = P.dimension();
    VT c_e(n);
    for(int i=0; i<n; i++){
        c_e(i)=res.first[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);

    Point xc = Point(n);
    //std::vector<Ball> S0;
    vecBalls.push_back(Ball(xc, res.second * res.second));
    ConvSet.push_back(Interballs(P.dimension(), vecBalls));
}*/

template <class VPolytope, class Interballs, class Ball,  typename FT, class Parameters>
void get_next_convex(VPolytope &VP, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls,
                     FT p_value, FT epsilon, Parameters &var, bool &done){

    typedef typename VPolytope::PolytopePoint Point;
    typedef typename Parameters::RNGType RNGType;
    Parameters var2 = var;
    var2.coordinate = true;
    int n = var.n;
    bool print = var.verbose;
    //if (print) std::cout<<"computation of "<<ConvSet.size()+1<<" encosing convex body started...\n"<<std::endl;
    Interballs Si = ConvSet[ConvSet.size()-1];
    //if (print) std::cout<<"last convex body has = "<<Si.num_of_balls()<<" balls"<<std::endl;
    std::list<Point> randPoints;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType &rng2 = var.rng;

    Ball test_ball;//= Si.get_ball();
    Point q(n);
    //if (print) std::cout<<"origin belongs to last convex body: "<<Si.is_in(q)<<std::endl;
    FT rad;
    if (ConvSet.size() ==1) {
        rad = vecBalls[0].radius();
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere<RNGType ,Point >(n, rad));
        }
    } else {
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var);
    }
    //q.print();
    //if (print) std::cout<<"q after belongs to last convex body: "<<Si.is_in(q)<<std::endl;

    std::list<Point> SirandPoints;
    int count = 0;
    int thre = int(3.0/(1.0-epsilon));
    bool too_few = false;

    while (true) {

        q = Point(n);
        do {

            //std::cout<<"count = "<<count<<" epsilon = "<<epsilon<<std::endl;
            if (count >= thre) {
                if (Si.num_of_balls() > ConvSet[ConvSet.size() - 1].num_of_balls()) {
                    q = Point(n);
                    SirandPoints.clear();
                    rand_point_generator(Si, q, 1200, var.walk_steps, SirandPoints, var2);
                    check_converg2<Point>(VP, SirandPoints, epsilon, done, too_few, var, false);
                } else {
                    check_converg2<Point>(VP, randPoints, epsilon, done, too_few, var, false);
                }
                //std::cout<<"VP convergence.... "<<"doneIter = "<<done<< " too_few = "<<too_few<<std::endl;
                if (done) {
                    if (Si.num_of_balls() > ConvSet[ConvSet.size() - 1].num_of_balls()) {
                        ConvSet.push_back(Si);
                    }
                    std::cout << "number of balls in Si = " << Si.num_of_hyperplanes() << std::endl;
                    return;

                } else if (!too_few) {
                    if (Si.num_of_balls() > ConvSet[ConvSet.size() - 1].num_of_balls()) {
                        ConvSet.push_back(Si);
                    }
                    std::cout << "number of balls in Si = " << Si.num_of_balls() << std::endl;
                    done = true;
                    return;
                } else {
                    count = -1000;
                    too_few = false;
                }
            }
            count++;
            rand_point(Si, q, var);
        } while (VP.is_in(q) == -1);

        construct_ball(Si, VP, ConvSet, randPoints, q, p_value, var);
        std::cout<<"number of balls = "<<Si.num_of_balls()<<std::endl;
        count = 0;
        too_few = false;
    }


}


template <class Polytope, class Interballs, class Ball, typename NT, class Parameters>
void get_ball_schedule(Polytope &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls,
                       NT p_value, NT epsilon, Parameters &var) {

    int n = var.n;
    typedef BallIntersectPolytope<Polytope, NT>        BallPoly;
    typedef typename Polytope::PolytopePoint Point;
    typedef typename Parameters::RNGType RNGType;
    bool print = var.verbose;
    //var.walk_steps=11;
    if (print) std::cout<<"computing first convex enclosing body...(minimum enclosing ball)\n"<<std::endl;
    get_first_conv(P, ConvSet, vecBalls, var);
    ConvSet.push_back(Interballs(P.dimension(), vecBalls));
    if (print) std::cout<<"P in first ball "<<P.is_in_ball(vecBalls[0].center(), vecBalls[0].radius())<<std::endl;
    P.print();
    if (print) std::cout<<"first convex enclosing body computed!\n"<<std::endl;
    if (print) std::cout<<"center = ";
    if (print) vecBalls[0].center().print();
    if (print) std::cout<<"radius = "<<vecBalls[0].radius()<<std::endl;
    Point p = get_point_in_Dsphere<RNGType , Point>(n, vecBalls[0].radius());
    if (print) std::cout<<"checking if first enclosing body is the last as well..\n"<<std::endl;
    //if (is_last_conv(P, ConvSet[0], a, var, last_ratio, true)) {
        //if (print) std::cout<<"firsti is last as well\n"<<std::endl;
        //return;
    //}
    //if (print) std::cout<<"first not last.... Computing sequence of enclosing convex bodies..\n"<<std::endl;
    bool done = false;
    int test_counter = 2;
    //while(!done) {
        get_next_convex(P, ConvSet, vecBalls, p_value, epsilon, var, done);
        //if (print) std::cout<<"computation of "<<ConvSet.size()<<" encosing convex body completed.."<<std::endl;
    //}
}

#endif
