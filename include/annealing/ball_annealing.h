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

/*
template <class T1, class T2, typename FT>
void reconstruct_ball(T1 &P, T2 &randPoints, Interballs &S, Interballs &Si, Point &center,
                      Point &q, Point &direction, FT a, std::vector<Interballs> &Convset, bool &last_ball, bool &last_conv, vars &var) {

    int n = var.n;
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    std::pair<NT,NT> bpair = Si.line_intersect(Point(n),v);
    FT min_plus = bpair.first, tol = 0.000001;
    //std::cout<<"min_plus = "<<min_plus<<std::endl;
    Point q2 = v * min_plus;

    //Point midq(n);
    //Point midpoint(n);
    //Point c0(n);
    //FT dis = std::sqrt((q-q2).squared_length());
    //Point temp = direction;// * (min_plus);
            //temp = temp * (1.0 / std::sqrt(temp.squared_length()));
    Interballs SiIter(n);
    typename  std::list<Point>::iterator pit;
    T2 listIter2;
    //int count = 0;
    Ball ball_iter;
    Point p(n);
    int count;


    //std::cout<<"closest center found, count = "<<count<<std::endl;
    //center.print();
    //std::cout<<"radius = "<<std::sqrt((center - q).squared_length())<<std::endl;
    //std::cout<<"P is in ball: "<<P.is_in_ball(center,  std::sqrt((center-q).squared_length()) )<<std::endl;
    FT rad1 = std::sqrt((center - q).squared_length());
    FT smallrad = rad1;
    FT rad2 = 2*std::sqrt((center - q2).squared_length());
    FT midrad;
    //rad2 = 2 * rad1;


    Point Inicenter = center;

    int iters=0;
    while(true) {
        iters++;
        ball_iter = Ball(Inicenter, rad2*rad2);
        count = 0;
        listIter2 = randPoints;
        pit = listIter2.begin();
        while(pit!=listIter2.end()) {
            if (ball_iter.is_in(*pit)  && Si.is_in(*pit)) {
                count++;
                p = (*pit);
                pit=listIter2.erase(pit);
            } else {
                pit++;
            }
        }
        if (count>200){
            break;
        } else {
            if (iters >10){
                SiIter = Si;
                SiIter.add_ball(ball_iter);
                if (is_last_ball(S, SiIter, p, 0.1, var)) {
                    Convset.push_back(SiIter);
                    last_ball = true;
                    return;
                }
                if (is_last_conv(P, Si, 0.1, var, false)){
                    Convset.push_back(Si);
                    last_conv = true;
                    return;
                }
                return;
            }
            std::cout<<"rad1 = "<<rad1<<" rad2 = "<<rad2<<" count = "<<count<<std::endl;
            rad2 = 10*rad2;
        }
    }
    p = Point(n);

    bool secondball = false;
    bool done = false, firstit = true;

    while (!done) {
        //firstit = true;
    while(true) {
        midrad = (rad1 + rad2) / 2.0;

        listIter2 = randPoints;
        pit = listIter2.begin();

        ball_iter = Ball(Inicenter, midrad*midrad);
        count = 0;
        while(pit!=listIter2.end()) {
            if (ball_iter.is_in(*pit)  && Si.is_in(*pit)) {
                count++;
                p = (*pit);
                pit=listIter2.erase(pit);
            } else {
                pit++;
            }
        }

        if (count<190){
           // firstit = false;
            std::cout<<"RECONSTRUCTION... count = "<<count<<std::endl;
            SiIter = Si;
            SiIter.add_ball(ball_iter);
            if(is_last_ball(S, SiIter, p, a, var)) {
                //Si.add_ball(ball_iter);
                Convset.push_back(SiIter);
                std::cout<<"\nADD CONVEX No. "<<Convset.size()<<std::endl;
                S = SiIter;
                p = Point(n);
                randPoints.clear();
                rand_point_generator(S, p, 1200, var.walk_steps, randPoints, var);
                if (midrad == smallrad) {
                    std::cout <<"MID==1ST TIGHT<=190"<<std::endl;
                    last_ball = true;
                    done = true;
                    break;
                }
                rad2 = midrad;
                rad1 = 2.0 * smallrad - rad2;
                //secondball = true;
                break;
            }
            rad1 = midrad;
            //q = midq;
        } else {
            if (midrad == smallrad) {
                //SiIter = Si;
                std::cout <<"MID==1ST TIGHT>190 count = "<<count<<std::endl;
                Si.add_ball(ball_iter);
                //Convset.push_back(Si);
                //std::cout<<"\nADD CONVEX No. "<<Convset.size()<<std::endl;
                //S = Si;
                //p = Point(n);
                //randPoints.clear();
                //rand_point_generator(S, p, 1200, var.walk_steps, randPoints, var);
                done = true;
                break;
            }
            std::cout << "No RECONSTRUCTION... count = " << count << std::endl;
            rad2 = midrad;
        }
        count = 0;

    }
    }



}*/


template <class Ball, class Polytope, class Point, class Parameters>
Ball construct_ball(Polytope &P, Point q, Point &cent, Point &qbound, Point &direction, Parameters &var) {

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
        }
    }
    VT a = Mat.colPivHouseholderQr().solve(b);
    if (a.dot(p) > 1.0) {
        a = -a;
        z0 = -1.0;
    }

    q = v * (min_plus + 0.0001);
    qbound = q;
    typename std::vector<NT>::iterator tempit = temp.iter_begin();
    int i = 0;
    for ( ; tempit!=temp.iter_end(); ++tempit, ++i) {
        *tempit = a(i);
    }
    direction = temp;
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
        if (P.is_in_ball(center, std::sqrt( (center - q).squared_length()) )) {
            break;
        }
    }
    std::cout<<"first ball found, counter = "<<counter<<std::endl;
    Point midpoint(n);
    counter = 0;
    std::cout<<"bisection method to find closest center\n";
    while ( std::sqrt((center - c0).squared_length())>tol ) {
        counter++;
        midpoint = (center + c0) * 0.5;
        if (P.is_in_ball(midpoint, std::sqrt( (midpoint - q).squared_length()) )) {
            center = midpoint;
        } else {
            c0 = midpoint;
        }
    }
    std::cout<<"closest center found, counter = "<<counter<<std::endl;
    center.print();
    std::cout<<"radius = "<<std::sqrt((center - q).squared_length())<<std::endl;
    cent = center;

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


template <class Interballs, class Point, typename FT, class Parameters>
bool is_last_ball(Interballs &S,Interballs &Si, Point &p, FT a, Parameters &var) {

    int n = var.n;
    bool print = var.verbose, check = false;
    //Interballs Si = (*ConvSet.end());
    std::list<Point> randPoints;

    rand_point_generator(S, p, 1200, var.walk_steps, randPoints, var);

    std::vector<FT> ratios;
    FT countsIn = 0.0;
    int i = 1;

    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (Si.is_in(*pit)) {
            countsIn += 1.0;
        }
        if (i % 120 == 0) {
            //if (print) std::cout<<"conv2conv ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countsIn/120.0);
            countsIn = 0.0;
        }
    }

    if (check_t_test(ratios, a, 0.1)) {
        return true;
    }
    return false;

}


template <class Polytope, class Interballs, typename FT, class Parameters>
bool is_last_conv(Polytope &P,Interballs &Si, FT a, Parameters &var, FT &last_ratio, bool first) {

    typedef typename Polytope::PolytopePoint Point;
    typedef typename Parameters::RNGType RNGType;
    int n = var.n;
    bool print = var.verbose, check = false;
    //Interballs Si = (*ConvSet.end());
    std::list<Point> randPoints;
    //Point p(n);
    Point p(n);

    FT rad;
    if (first) {
        if (print) std::cout<<"this is first ball\n"<<std::endl;
        rad = Si.get_ball(0).radius();
        for (int i = 0; i < 10; ++i) {
            randPoints.push_back(get_point_in_Dsphere<RNGType, Point >(n, rad));
        }
    } else {
        rand_point_generator(Si, p, 10, var.walk_steps, randPoints, var);
    }

    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
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
            randPoints.push_back(get_point_in_Dsphere<RNGType ,Point >(n, rad));
        }
    }else {
        rand_point_generator(Si, p, 1200, var.walk_steps, randPoints, var);
    }
    std::vector<FT> ratios;
    FT countsIn = 0.0;
    int i = 1;

    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (P.is_in(*pit)) {
            countsIn += 1.0;
        }
        if (i % 120 == 0) {
            if (print) std::cout<<"ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countsIn/120.0);
            countsIn = 0.0;
        }
    }

    std::pair<FT,FT> mv;
    if (check_t_test(ratios, a, 0.1)) {
        mv = getMeanVariance(ratios);
        last_ratio = mv.first;
        return true;
    }
    mv = getMeanVariance(ratios);
    last_ratio = mv.first;
    return false;

}

template <class Polytope, class PointList, class Interballs, class Point, typename NT, class Parameters>
void reconstruct_ball(Polytope &P, PointList &randPoints, Interballs &S, Interballs &Si, Point &center,
                      Point &q, Point &direction, NT a, std::vector<Interballs> &Convset, bool &last_ball,
                      bool &last_conv, NT &last_ratio, Parameters &var) {

    int n = var.n;
    typedef typename Interballs::cball Ball;
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    std::pair<NT,NT> bpair = Si.line_intersect(Point(n),v);
    NT min_plus = bpair.first, tol = 0.000001;
    //std::cout<<"min_plus = "<<min_plus<<std::endl;
    Point q2 = v * min_plus;

    //Point midq(n);
    //Point midpoint(n);
    //Point c0(n);
    //FT dis = std::sqrt((q-q2).squared_length());
    //Point temp = direction;// * (min_plus);
    //temp = temp * (1.0 / std::sqrt(temp.squared_length()));
    Interballs SiIter(n);
    typename  std::list<Point>::iterator pit;
    PointList listIter2;
    //int count = 0;
    Ball ball_iter;
    Point p(n);
    int count;


    //std::cout<<"closest center found, count = "<<count<<std::endl;
    //center.print();
    //std::cout<<"radius = "<<std::sqrt((center - q).squared_length())<<std::endl;
    //std::cout<<"P is in ball: "<<P.is_in_ball(center,  std::sqrt((center-q).squared_length()) )<<std::endl;
    NT rad1 = std::sqrt((center - q).squared_length());
    NT smallrad = rad1;
    //FT rad2 = 2*std::sqrt((center - q2).squared_length());
    NT midrad;
    NT rad2 = 2 * rad1;


    Point Inicenter = center;

    int iters=0;
    while(true) {
        iters++;
        ball_iter = Ball(Inicenter, rad2*rad2);
        count = 0;
        listIter2 = randPoints;
        pit = listIter2.begin();
        while(pit!=listIter2.end()) {
            if (ball_iter.is_in(*pit)  && Si.is_in(*pit)) {
                count++;
                p = (*pit);
                pit=listIter2.erase(pit);
            } else {
                pit++;
            }
        }
        if (count>190){
            break;
        } else {
            if (iters >10){
                SiIter = Si;
                SiIter.add_ball(ball_iter);
                if (is_last_ball(S, SiIter, p, 0.1, var)) {
                    Convset.push_back(SiIter);
                    last_ball = true;
                    return;
                }
                if (is_last_conv(P, Si, 0.1, var, last_ratio, false)){
                    Convset.push_back(Si);
                    last_conv = true;
                    return;
                }
                Convset.push_back(SiIter);
                last_ball = true;
                return;
                //return;
            }
            std::cout<<"rad1 = "<<rad1<<" rad2 = "<<rad2<<" count = "<<count<<std::endl;
            rad2 = 10*rad2;
        }
    }
    p = Point(n);

    bool secondball = false;
    bool done = false, firstit = true;

    while (!done) {
        //firstit = true;
        while(true) {
            midrad = (rad1 + rad2) / 2.0;

            listIter2 = randPoints;
            pit = listIter2.begin();

            ball_iter = Ball(Inicenter, midrad*midrad);
            count = 0;
            while(pit!=listIter2.end()) {
                if (ball_iter.is_in(*pit)  && Si.is_in(*pit)) {
                    count++;
                    p = (*pit);
                    pit=listIter2.erase(pit);
                } else {
                    pit++;
                }
            }

            if (count<190){
                // firstit = false;
                std::cout<<"RECONSTRUCTION... count = "<<count<<std::endl;
                SiIter = Si;
                SiIter.add_ball(ball_iter);
                if(is_last_ball(S, SiIter, p, a, var)) {
                    //Si.add_ball(ball_iter);
                    Convset.push_back(SiIter);
                    std::cout<<"\nADD CONVEX No. "<<Convset.size()<<std::endl;
                    S = SiIter;
                    p = Point(n);
                    randPoints.clear();
                    rand_point_generator(S, p, 1200, var.walk_steps, randPoints, var);
                    if (midrad == smallrad) {
                        std::cout <<"MID==1ST TIGHT<=190"<<std::endl;
                        last_ball = true;
                        done = true;
                        break;
                    }
                    rad2 = midrad;
                    rad1 = 2.0 * smallrad - rad2;
                    //secondball = true;
                    break;
                }
                rad1 = midrad;
                //q = midq;
            } else {
                if (midrad == smallrad) {
                    //SiIter = Si;
                    std::cout <<"MID==1ST TIGHT>190 count = "<<count<<std::endl;
                    Si.add_ball(ball_iter);
                    //Convset.push_back(Si);
                    //std::cout<<"\nADD CONVEX No. "<<Convset.size()<<std::endl;
                    //S = Si;
                    //p = Point(n);
                    //randPoints.clear();
                    //rand_point_generator(S, p, 1200, var.walk_steps, randPoints, var);
                    done = true;
                    break;
                }
                std::cout << "No RECONSTRUCTION... count = " << count << std::endl;
                rad2 = midrad;
            }
            count = 0;

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
    std::cout<<"radius of minim ball = "<<rad<<std::endl;
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

template <class Polytope, class Interballs, class Ball,  typename FT, class Parameters>
void get_next_convex(Polytope &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls, FT a, Parameters &var, FT &last_ratio, bool &done){

    typedef typename Polytope::PolytopePoint Point;
    typedef typename Parameters::RNGType RNGType;
    int n = var.n;
    bool print = var.verbose;
    if (print) std::cout<<"computation of "<<ConvSet.size()+1<<" encosing convex body started...\n"<<std::endl;
    Interballs Si = ConvSSi, VP, PolyBallSet, randPoints, added_in_set, varet[ConvSet.size()-1];
    if (print) std::cout<<"last convex body has = "<<Si.num_of_balls()<<" balls"<<std::endl;
    std::list<Point> randPoints;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType &rng2 = var.rng;

    Ball test_ball;//= Si.get_ball();
    for (int j = 0; j < Si.num_of_balls(); ++j) {

    }
    Point q(n);
    if (print) std::cout<<"origin belongs to last convex body: "<<Si.is_in(q)<<std::endl;
    FT rad;
    if (ConvSet.size() ==1) {
        rad = vecBalls[0].radius();
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere<RNGType ,Point >(n, rad));
        }
    } else {
        rand_point_generator(Si, q, 1200, var.walk_steps, randPoints, var);
    }
    q.print();
    if (print) std::cout<<"q after belongs to last convex body: "<<Si.is_in(q)<<std::endl;

    std::list<Point> listIter = randPoints;
    std::list<Point> listIter2;

    //int rand_coord;
    Ball ball_iter;
    //std::vector<int> index(1200, -1);
    int count = 0, count_bef=1200;
    FT totcount = 0, itercount = 0;
    Interballs SiIter(n);
    Interballs S = Si;
    Point direction(n), cent(n), qbound(n);

    bool added=false;
    bool tested=false, last_ball = false, last_conv = false;
    typename  std::list<Point>::iterator pit;
    int too_few=0;
        while (true) {

            //count_bef = count;
            if (print) std::cout<<"points outside new convex = "<<listIter.size()<<std::endl;
            //rand_coord = uidist(rng);
            //index.assign(1200,-1);
            if (print) std::cout<<"q_bef  is in: "<<Si.is_in(q)<<std::endl;
            if (ConvSet.size() ==1 && Si.num_of_balls()==1) {
                q = get_point_in_Dsphere<RNGType ,Point >(n, rad);
            } else {
                rand_point(Si, q, var);
            }
            //totcount += 0.0;
            if (print) std::cout<<"rand point is in: "<<Si.is_in(q)<<std::endl;
            //q.print();
            /*if (is_last_conv(P, Si, a, var, false)){
                ConvSet.push_back(Si);
                done = true;
                return;
            }*/
            if (!P.is_in(q)) {
                totcount=0.0;
                if (print) std::cout<<"construct ball"<<std::endl;
                ball_iter = construct_ball<Ball>(P, q, cent, qbound, direction, var);
                if (print) std::cout<<"ball constructed"<<std::endl;
                listIter2 = randPoints;
                pit = listIter2.begin();
                while(pit!=listIter2.end()) {
                    if (ball_iter.is_in(*pit)  && Si.is_in(*pit)) {
                        count++;
                        pit=listIter2.erase(pit);
                    } else {
                        pit++;
                    }
                }
                listIter = listIter2;
                vecBalls.push_back(ball_iter);


                itercount = 0.0;
                //totcount = 0.0;
                std::cout<<"count = "<<count<<std::endl;
                if (count<190) {
                    std::cout<<"count = "<<count<<std::endl;
                    SiIter = Si;
                    SiIter.add_ball(ball_iter);
                    if (is_last_ball(S, SiIter, q, a, var)) {
                        Si.add_ball(ball_iter);
                        ConvSet.push_back(Si);
                        std::cout<<"\nADD CONVEX No. "<<ConvSet.size()<<std::endl;
                        return;
                    }
                    reconstruct_ball(P, randPoints, S, Si, cent, qbound, direction, a,
                                     ConvSet, last_ball, last_conv, last_ratio, var);
                    if (last_conv){
                        done = true;
                        return;
                    }
                    if (last_ball) return;
                    //ConvSet.push_back(Si);
                    //return;
                    q = Point(n);
                    totcount=0.0;
                    continue;
                } else if (count_bef - count < 10 && false) {
                    std::cout<<"count_bef - count = "<<count_bef - count<<std::endl;
                    if (is_last_conv(P, Si, a, var, last_ratio, false)){
                        ConvSet.push_back(Si);
                        done = true;
                        return;
                    }
                }

                if (added) {
                    if (count_bef - count < 12) {
                        too_few++;
                        std::cout << "too_few = " << too_few << std::endl;
                        std::cout << "count_bef - count = " << count_bef - count << std::endl;

                        if (too_few >= 150) {
                            //continue;
                            std::cout << "count_bef - count = " << count_bef - count << std::endl;
                            //if (is_last_conv(P, Si, a, var, last_ratio, false) || true) {
                                ConvSet.push_back(Si);
                                std::cout << "\nADD LAST CONVEX No. " << ConvSet.size() << std::endl;
                                done = true;
                                return;
                            //} else {
                                added = false;
                            //}
                        }
                        if (count_bef - count == 0) {
                            count=0;
                            totcount=0.0;
                            continue;
                        }
                    } else {
                        too_few = 0;
                    }
                }
                q = Point(n);
                Si.add_ball(ball_iter);
                added = true;
                tested = false;
                totcount =0.0;
                count_bef = count;
                count = 0;
                std::cout<<"BAL ADDED No: "<<Si.num_of_balls()<<std::endl;
            } else {
                itercount += 1.0;
                totcount += 1.0;
                if (totcount>3.0 && !tested) {
                    std::cout<<"totcount = "<<totcount<<std::endl;
                    if (is_last_conv(P, Si, a, var, last_ratio, false)){
                        ConvSet.push_back(Si);
                        done = true;
                        return;
                    }
                    tested = true;
                }
            }
        }
    ConvSet.push_back(Si);
}


template <class Polytope, class Interballs, class Ball, typename NT, class Parameters>
void get_ball_schedule(Polytope &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls,
                       NT a, NT &last_ratio, Parameters &var) {

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
    if (is_last_conv(P, ConvSet[0], a, var, last_ratio, true)) {
        if (print) std::cout<<"firsti is last as well\n"<<std::endl;
        return;
    }
    if (print) std::cout<<"first not last.... Computing sequence of enclosing convex bodies..\n"<<std::endl;
    bool done = false;
    int test_counter = 2;
    while(!done) {
        get_next_convex(P, ConvSet, vecBalls, a, var, last_ratio, done);
        if (print) std::cout<<"computation of "<<ConvSet.size()<<" encosing convex body completed.."<<std::endl;
    }
}
