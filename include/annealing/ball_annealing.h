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
    NT min_plus = intersect_line_Vpoly(P.get_mat(), center, v, false);
    q = v * (min_plus * 1.01);
    //bool isin = false;
    Point c0(n);
    NT const tol = 0.000001;

    // find the first ball that contains v-polytope P
    while (true) {
        c0 = center;
        center = center - (min_plus * v);
        if (P.is_in_ball(center, std::sqrt( (center - q).squared_length()) )) {
            break;
        }
    }
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

    if (p_mval > p_test + t_value*(ni-1)*(p_varval/std::sqrt(FT(ni)))) {
        passed = true;
    }
    return passed;
}


template <class T, typename FT>
bool is_last_conv(T &P, std::vector<Interballs> &ConvSet, Point &p, FT a, vars &var) {

    int n = var.n;
    bool print = var.verbose, check = false;
    Interballs Si = (*ConvSet.end());
    std::list<Point> randPoints;
    //Point p(n);

    rand_point_generator(Si, p, 10, var.walk_steps, randPoints, var);

    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        if (P.is_in(*pit)) {
            check = true;
            break;
        }
    }
    if (!check) return false;

    randPoints.clear();
    rand_point_generator(Si, p, 1200, var.walk_steps, randPoints, var);
    std::vector<FT> ratios;
    FT countsIn = 0.0;
    int i = 1;

    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (P.is_in(*pit)) {
            countIn += 1.0;
        }
        if (i % 120 == 0) {
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
void get_first_conv(T &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls) {
    std::pair<Point, NT> res = compute_minball(P);
    int n = P.dimension();
    VT c_e(n);
    for(int i=0; i<n; i++){
        c_e(i)=res.first[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);

    xc = Point(n);
    //std::vector<Ball> S0;
    vecBalls.push_back(Ball(xc, res.second * res.second));
    ConvSet.push_back(InterBalls(P.dimension(), vecBalls));
}

template <class T1, typename FT>
void get_next_convex(T1 &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls, Point &p, FT a, vars &var, bool &check){

    int n =var.n;
    Interballs Si = (*ConvSet.end());
    std::list<Point> randPoints;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType &rng2 = var.rng;

    if (ConvSet.size() ==1) {
        FT rad = vecBalls[0].radius();
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere(n, rad));
        }
    } else {
        rand_point_generator(Si, p, 1200, var.walk_steps, randPoints, var);
    }

    std::list<Point> listIter = randPoints;
    std::list<Point> listIter2;
    Point q(n);// = rand_point(Si, p, )
    int rand_coord;
    Ball ball_iter;
    //std::vector<int> index(1200, -1);
    int count = 0;
    typename  std::list<Point>::iterator pit;
    while (listIter.size() > 190) {
        //rand_coord = uidist(rng);
        //index.assign(1200,-1);
        rand_point(Si, q, var);
        if (!P.is_in(q)) {
            ball_iter = construct_ball(P, q, var);
            listIter2 = listIter;
            pit = listIter2.begin();
            while(pit!=listIter2.end()) {
                if (!ball_iter.is_in(*pit)) {
                    count++;
                    pit=listIter2.erase(pit);
                } else {
                    pit++;
                }
            }
            if (count>=60) {
                listIter = listIter2;
                vecBalls.push_back(ball_iter);
                Si.add_ball(ball_iter);
                q = Point(n);
            }
            count = 0;
        }
    }
    ConvSet.push_back(Si);

}


template <class T1, typename FT>
void get_ball_schedule(T1 &P, std::vector<Interballs> &ConvSet, std::vector<Ball> &vecBalls, FT a, vars &var) {

    int n = var.n;
    get_first_conv(P, ConvSet, vecBalls);
    Point p = get_point_in_Dsphere(n, vecBalls[0].radius());
    if (is_last_conv(P, ConvSet, p, a, var)) {
        return;
    }
    bool check = false;
    while(true) {
        get_next_convex(P, ConvSet, vecBalls, p, a, var, check);
        if (check) {
            if (is_last_conv(P, p, a, var)) return;
        }
    }
}