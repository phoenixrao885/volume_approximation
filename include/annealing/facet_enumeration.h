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

#ifndef FACET_ENUMERATION_H
#define FACET_ENUMERATION_H

template <class Polytope, class Ball, class Point, class Parameters>
void enclosing_ball(Polytope &P, Ball &B0, Point &center, Parameters &var) {

    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Polytope::NT NT;
    //typedef typename Polytope::PolytopePoint Point;
    unsigned int n = P.dimension();
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
    std::cout<<"radius of minim ball = "<<rad<<std::endl;
    xc.print();
    center = Point(n) - xc;
    //Point xc = Point(n);
    //std::vector<Ball> S0;
    B0 = Ball(Point(n), rad*rad);
    //ConvSet.push_back(Interballs(P.dimension(), vecBalls));

}


template <class HPolytope, class VPolytope, typename NT, class Parameters>
void check_convergence(HPolytope &HP, VPolytope &VP, NT epsilon, bool &done, Parameters &var){

    typedef typename VPolytope::PolytopePoint Point;


    int M = 1200;
    int ni = 10;
    int Mni = 120;
    int n = var.n;
    NT a = 0.05;
    std::list<Point> randPoints;
    Point p(n);

    rand_point_generator(HP, p, M, var.walk_steps, randPoints, var);

    std::vector<NT> ratios;
    NT countsIn = 0.0;
    int i = 1;

    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, ++i){
        if (VP.is_in(*pit)) {
            countsIn += 1.0;
        }
        if (i % 120 == 0) {
            std::cout<<"conv2conv ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countsIn/120.0);
            countsIn = 0.0;
        }
    }

    if (check_t_test(ratios, a, epsilon)) {
        done = true;
    }

}


template <class ConvexBody, class VPolytope, class Point>
void add_facet(ConvexBody &K, VPolytope &VP, Point &q){

    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::VT VT;

    unsigned int n = VP.dimension(), i;
    NT z0;
    VT a(n);
    Point center(n);
    std::vector<NT> hyp(n,0);
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    NT min_plus = intersect_line_Vpoly<NT>(VP.get_mat(), center, v, false, false);
    q = v * (min_plus + 0.0001);
    get_separeting_hyp(VP.get_mat(), q, hyp, z0);

    typename std::vector<NT>::iterator pit;
    pit = hyp.begin(); i= 0;
    for ( ; pit!=hyp.end(); ++pit, ++i) {
        a(i) = *pit;
    }
    K.add_facet(a, z0);
}

template <class ConvexBody, class HPolytope, class VPolytope, class Point>
void add_facet(ConvexBody &K, HPolytope &HP, VPolytope &VP, Point &q){

    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::VT VT;

    unsigned int n = VP.dimension(), i;
    NT z0;
    VT a(n);
    Point center(n);
    std::vector<NT> hyp(n,0);
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    NT min_plus = intersect_line_Vpoly<NT>(VP.get_mat(), center, v, false, false);
    q = v * (min_plus + 0.0001);
    get_separeting_hyp(VP.get_mat(), q, hyp, z0);

    typename std::vector<NT>::iterator pit;
    pit = hyp.begin(); i= 0;
    for ( ; pit!=hyp.end(); ++pit, ++i) {
        a(i) = *pit;
    }
    K.add_facet(a, z0);
    HP.add_facet(a,z0);
}


template <class ConvexBody, class VPolytope, class Point>
void add_facet2(ConvexBody &K, VPolytope &VP, int & on_faces, Point &q){

    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::VT VT;
    typedef typename VPolytope::MT MT;

    unsigned int n = VP.dimension(), i;
    std::vector <NT> lambdas(VP.num_of_vertices());
    NT z0 = 1.0;
    //VT a(n);
    Point center(n);
    //std::vector <NT> hyp(n, 0);
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    NT min_plus = intersect_line_Vpoly2<NT>(VP.get_mat(), center, v, false, lambdas);
    MT Mat = MT::Zero(n,n);
    VT b = VT::Ones(n);
    VT p(n);
    MT V = VP.get_mat();
    int count = 0;

    for (int j = 0; j < VP.num_of_vertices(); ++j) {
        if (lambdas[j] > 0.0) {
            Mat.row(count) = V.row(j);
            count++;
        } else {
            p = V.row(j);
        }
    }
    if (count < n) {
        on_faces++;
        return;
    }
    VT a = Mat.colPivHouseholderQr().solve(b);
    if (a.dot(p) > 1.0) {
        a = -a;
        z0 = -1.0;
    }

    K.add_facet(a, z0);
}

template <class ConvexBody, class HPolytope, class VPolytope, class Point>
void add_facet2(ConvexBody &K, HPolytope &HP, VPolytope &VP, bool &added, Point &q) {

    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::VT VT;
    typedef typename VPolytope::MT MT;

    unsigned int n = VP.dimension(), i;
    std::vector <NT> lambdas(VP.num_of_vertices());
    NT z0 = 1.0;
    //VT a(n);
    Point center(n);
    //std::vector <NT> hyp(n, 0);
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    NT min_plus = intersect_line_Vpoly2<NT>(VP.get_mat(), center, v, false, lambdas);
    MT Mat = MT::Zero(n,n);
    VT b = VT::Ones(n);
    VT p(n);
    MT V = VP.get_mat();
    int count = 0;

    for (int j = 0; j < VP.num_of_vertices(); ++j) {
        if (lambdas[j] > 0.0) {
            Mat.row(count) = V.row(j);
            count++;
        } else {
            p = V.row(j);
        }
    }
    if(count<n){
        //std::cout<<count<<std::endl;
        return;
    }
    added = true;
    VT a = Mat.colPivHouseholderQr().solve(b);
    if (a.dot(p) > 1.0) {
        a = -a;
        z0 = -1.0;
    }

    K.add_facet(a, z0);
    HP.add_facet(a, z0);
}


template <class ConvexBody, class HPolytope, class VPolytope, typename NT, class Parameters>
void construct_simplex(ConvexBody &K, HPolytope &HP, VPolytope &VP, NT epsilon, bool &done, int &on_faces, Parameters &var){

    typedef typename ConvexBody::VT VT;
    //typedef typename VPolytope::NT NT;
    typedef typename VPolytope::PolytopePoint Point;
    typedef typename Parameters::RNGType RNGType;

    unsigned int n = K.dimension(), i;
    //std::pair<VT, NT> halfspace;
    Point q(n);
    int count, j = 0;
    int thre = int(4.0/(1.0-epsilon));
    bool added = false;

    while(j<2*n*n) {
    //for (int j = 0; j < 2*n; ++j) {
        q = Point(n);
        count = 0;

        do{
            if (count>=thre) {
                check_convergence(HP, VP, epsilon, done, var);
                if (done){
                    return;
                } else {
                    count=-1000;
                }
            }
            count++;
            //std::cout<<K.num_of_hyperplanes()<<std::endl;
            if (K.num_of_hyperplanes() == 0) {
                q = get_point_in_Dsphere<RNGType, Point>(n, K.second().radius());
            } else {
                rand_point(K, q, var);
            }
        } while (VP.is_in(q)==-1);

        add_facet2(K, HP, VP, added, q);
        if (added) {
            std::cout<<added<<std::endl;
            added = false;
            j++;
        } else {
            on_faces++;
        }
        //add_facet(K, HP, VP, q);
        //j++;

    }

};

template <class HPolytope, class VPolytope, typename NT, class Parameters>
void detect_facet(HPolytope &HP, VPolytope &VP, NT epsilon, bool &done, int &on_faces, Parameters &var){

    typedef typename HPolytope::VT VT;
    typedef typename VPolytope::PolytopePoint Point;

    //NT z0;
    unsigned int n = HP.dimension(), i;
    //std::pair<VT, NT> halfspace;
    Point q(n);
    //VT a(n);
    //Point center(n);
    //std::vector<NT> hyp(n,0);
    int count = 0;
    int thre = int(4.0/(1.0-epsilon));

    do{
        if (count>=thre) {
            check_convergence(HP, VP, epsilon, done, var);
            if (done){
                return;
            } else {
                count=-1000;
            }
        }
        count++;
        rand_point(HP, q, var);
    } while (VP.is_in(q)==-1);

    add_facet2(HP, VP, on_faces, q);
    std::cout<<"number of facets = "<<HP.num_of_hyperplanes()<<std::endl;

};


template <class HPolytope, class Ball, class VPolytope, typename NT, class Parameters>
HPolytope facet_enumeration(VPolytope &VP, NT &epsilon, Parameters &var) {

    typedef BallIntersectPolytope<HPolytope,Ball> BallPoly;
    typedef typename Parameters::RNGType RNGType;
    typedef typename VPolytope::PolytopePoint Point;
    typedef typename VPolytope::VT VT;

    unsigned int n = var.n;
    HPolytope HP(n);
    Point center(n);
    bool done = false, check = false;
    int on_faces = 0;
    Ball B0;
    enclosing_ball(VP, B0, center, var);
    BallPoly BP(HP, B0);
    construct_simplex(BP, HP, VP, epsilon, done, on_faces, var);
    if(done) {
        done = false;
    }

    std::cout<<"number of facets after simplex = "<<HP.num_of_hyperplanes()<<std::endl;
    while (!done) {
        detect_facet(HP, VP, epsilon, done, on_faces, var);
    }

    //std::cout<<"nuber of facets after simplex = "<<BP.num_of_hyperplanes()<<std::endl;

    center.print();
    VT c_e(n);
    for (int i = 0; i < n; ++i) {
        c_e(i) = center[i];
    }
    HP.shift(c_e);
    HP.normalization();
    std::cout<<"on faces = "<<on_faces<<std::endl;
    HP.print();
    return HP;

};

#endif
