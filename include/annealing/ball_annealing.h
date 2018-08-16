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

template <class T>
Ball construct_ball(T &P, Point q, vars &var) {

    int n = var.n;
    Point center(n);
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    std::cout<<"v = ";
    v.print();
    NT min_plus = intersect_line_Vpoly(P.get_mat(), center, v, false);
    std::cout<<"min_plus = "<<min_plus<<std::endl;
    q.print();
    q = v * (min_plus * 1.01);
    std::cout<<"q = ";
    q.print();
    //std::cout<<"q is in P = "<<P.is_in(q);
    //bool isin = false;
    Point c0(n);
    NT const tol = 0.000001;

    // find the first ball that contains v-polytope P
    std::cout<<"find first ball\n";
    Point temp(n);
    temp = q;
    while (true) {
        c0 = center;

        center = center - q;
        std::cout<<std::sqrt( (center - q).squared_length())<<std::endl;
        if (P.is_in_ball(center, std::sqrt( (center - q).squared_length()) )) {
            break;
        }
    }
    std::cout<<"first ball found"<<std::endl;
    Point midpoint(n);
    while ( std::sqrt((center - c0).squared_length())>tol ) {
        midpoint = (center + c0) * 0.5;
        if (P.is_in_ball(midpoint, std::sqrt( (midpoint - q).squared_length()) )) {
            center = midpoint;
        } else {
            c0 = midpoint;
        }
    }

    return Ball(center,  (center - q).squared_length() );
}

template <typename FT>
bool check_t_test(std::vector<FT> ratios, FT a, FT p_test){

    bool passed = false;
    std::pair<FT,FT> mv = getMeanVariance(ratios);
    FT t_value = 0.700;
    FT p_mval = mv.first;
    FT p_varval = mv.second;
    int ni = ratios.size();
    std::cout<<"mean = "<<mv.first<<" var = "<<mv.second<<std::endl;
    if(mv.first >=190.0/1200.0 && false) {
        return true;
    }

    std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(FT(ni)))<<std::endl;
    if (p_mval > p_test + t_value*(ni-1)*(p_varval/std::sqrt(FT(ni)))) {
        passed = true;
    }
    return passed;
    //return false;
}


template <class T, typename FT>
bool is_last_conv(T &P,Interballs &Si, Point &p, FT a, vars &var, bool first) {

    int n = var.n;
    bool print = var.verbose, check = false;
    //Interballs Si = (*ConvSet.end());
    std::list<Point> randPoints;
    //Point p(n);

    FT rad;
    if (first) {
        if (print) std::cout<<"this is first ball\n"<<std::endl;
        rad = Si.get_ball(0).radius();
        for (int i = 0; i < 10; ++i) {
            randPoints.push_back(get_point_in_Dsphere(n, rad));
        }
    } else {
        rand_point_generator(Si, p, 10, var.walk_steps, randPoints, var);
    }

    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        if (P.is_in(*pit)) {
            check = true;
            break;
        }
    }
    if (print) std::cout<<"check is = "<<check<<std::endl;
    if (!check) return false;

    randPoints.clear();
    //FT rad = vecBalls[0].radius();
    if (first) {
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere(n, rad));
        }
    }else {
        rand_point_generator(Si, p, 1200, var.walk_steps, randPoints, var);
    }
    std::vector<FT> ratios;
    FT countsIn = 0.0;
    int i = 1;

    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (P.is_in(*pit)) {
            countsIn += 1.0;
        }
        if (i % 120 == 0) {
            if (print) std::cout<<"ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countsIn/120.0);
            countsIn = 0.0;
        }
    }

    if (check_t_test(ratios, a, 0.1)) {
        return true;
    }
    return false;

}


template <class T>
void get_first_conv(T &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls, vars &var) {

    int n = P.dimension();
    std::vector<NT> vec(n,0.0);

    Point xc(n);

    MT V = P.get_mat();
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

    //Point xc = Point(n);
    //std::vector<Ball> S0;
    vecBalls.push_back(Ball(Point(n), rad));
    ConvSet.push_back(Interballs(P.dimension(), vecBalls));

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

template <class T1, typename FT>
void get_next_convex(T1 &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls, FT a, vars &var, bool &done){

    int n = var.n;
    bool print = var.verbose;
    if (print) std::cout<<"computation of "<<ConvSet.size()+1<<" encosing convex body started...\n"<<std::endl;
    Interballs Si = ConvSet[ConvSet.size()-1];
    std::list<Point> randPoints;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType &rng2 = var.rng;

    Point q(n);
    FT rad;
    if (ConvSet.size() ==1) {
        rad = vecBalls[0].radius();
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere(n, rad));
        }
    } else {
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var);
    }

    std::list<Point> listIter = randPoints;
    std::list<Point> listIter2;

    //int rand_coord;
    Ball ball_iter;
    //std::vector<int> index(1200, -1);
    int count = 0;
    FT totcount = 0, itercount = 0;
    typename  std::list<Point>::iterator pit;
        while (listIter.size() > 1010) {
            if (print) std::cout<<"points outside new convex = "<<listIter.size()<<std::endl;
            //rand_coord = uidist(rng);
            //index.assign(1200,-1);
            if (ConvSet.size() ==1) {
                q = get_point_in_Dsphere(n, rad);
            } else {
                rand_point(Si, q, var);
            }
            if (print) std::cout<<"rand point";
            q.print();
            if (!P.is_in(q)) {
                if (print) std::cout<<"construct ball"<<std::endl;
                ball_iter = construct_ball(P, q, var);
                if (print) std::cout<<"ball constructed"<<std::endl;
                listIter2 = listIter;
                pit = listIter2.begin();
                while(pit!=listIter2.end()) {
                    if (ball_iter.is_in(*pit)) {
                        count++;
                        pit=listIter2.erase(pit);
                    } else {
                        pit++;
                    }
                }
                if (count>=10) {
                    listIter = listIter2;
                    vecBalls.push_back(ball_iter);
                    Si.add_ball(ball_iter);
                    q = Point(n);
                    itercount = 0.0;
                    totcount = 0.0;
                } else {
                    totcount += 1.0;
                }
                count = 0;
            } else {
                itercount += 1.0;
                totcount += 1.0;
                if (itercount / totcount >0.1 && totcount>10.0) {
                    if (is_last_conv(P, Si, q, a, var, false)){
                        ConvSet.push_back(Si);
                        done = true;
                        return;
                    }
                }
            }
        }
    ConvSet.push_back(Si);
}


template <class T1, typename FT>
void get_ball_schedule(T1 &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls, FT a, FT &round_value, vars &var) {

    int n = var.n;
    bool print = var.verbose;
    if (print) std::cout<<"computing first convex enclosing body...(minimum enclosing ball)\n"<<std::endl;
    get_first_conv(P, ConvSet, vecBalls, var);
    if (print) std::cout<<"P in first ball "<<P.is_in_ball(vecBalls[0].center(), vecBalls[0].radius())<<std::endl;
    P.print();
    if (print) std::cout<<"first convex enclosing body computed!\n"<<std::endl;
    if (print) std::cout<<"center = ";
    if (print) vecBalls[0].center().print();
    if (print) std::cout<<"radius = "<<vecBalls[0].radius()<<std::endl;
    Point p = get_point_in_Dsphere(n, vecBalls[0].radius());
    if (print) std::cout<<"checking if first enclosing body is the last as well..\n"<<std::endl;
    if (is_last_conv(P, ConvSet[0], p, a, var, true)) {
        if (print) std::cout<<"firsti is last as well\n"<<std::endl;
        return;
    }
    if (print) std::cout<<"first not last.... Computing sequence of enclosing convex bodies..\n"<<std::endl;
    bool done = false;
    int test_counter = 2;
    while(!done) {
        get_next_convex(P, ConvSet, vecBalls, a, var, done);
        if (print) std::cout<<"computation of "<<ConvSet.size()<<" encosing convex body completed.."<<std::endl;
    }
}