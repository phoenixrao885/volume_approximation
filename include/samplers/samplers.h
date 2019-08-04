// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H


// Pick a random direction as a normilized vector
template <class RNGType, class Point, typename NT>
Point get_direction(unsigned int dim) {

    boost::normal_distribution<> rdist(0,1);
    std::vector<NT> Xs(dim,0);
    NT normal = NT(0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType rng2 = var.rng;
    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = rdist(rng);
        normal += Xs[i] * Xs[i];
    }
    normal=1.0/std::sqrt(normal);

    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = Xs[i] * normal;
    }
    Point p(dim, Xs.begin(), Xs.end());
    return p;
}


// Pick a random point from a d-sphere
template <class RNGType, class Point, typename NT>
Point get_point_on_Dsphere(unsigned int dim, NT radius){
    Point p = get_direction<RNGType, Point, NT>(dim);
    p = (radius == 1.0) ? p : radius * p;
    return p;
}


// Pick a random point from a d-ball
template <class RNGType, class Point, typename NT>
Point get_point_in_Dsphere(unsigned int dim, NT radius){

    boost::random::uniform_real_distribution<> urdist(0,1);
    NT U;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng2(seed);
    Point p = get_direction<RNGType, Point, NT>(dim);
    U = urdist(rng2);
    U = std::pow(U, 1.0/(NT(dim)));
    p = (radius*U)*p;
    return p;
}

// ball walk with uniform target distribution
template <class RNGType, class Point, class Polytope, typename NT>
void ball_walk(Point &p,
               Polytope &P,
               NT delta)
{
    //typedef typename Parameters::RNGType RNGType;
    Point y = get_point_in_Dsphere<RNGType, Point>(p.dimension(), delta);
    y = y + p;
    if (P.is_in(y)==-1) p = y;
}

// WARNING: USE ONLY WITH BIRKHOFF POLYOPES
// Compute more random points using symmetries of birkhoff polytope
/*
template <class T, class Point, class K, typename  NT>
int birk_sym(T &P, K &randPoints, Point &p) {
    int n=std::sqrt(p.dimension());
    std::vector<int> myints;
    for (unsigned int i=0; i<n; i++){
        myints.push_back(i);
    }

    //The n! possible permutations with n elements
    do {
        std::vector<NT> newpv;

        for (unsigned int j=0; j<p.dimension(); j++){
            int idx = (myints[j/n])*n+1+j%n-1;
            newpv.push_back(p[idx]);
        }

        Point new_p(p.dimension(),newpv.begin(),newpv.end());
        if(P.is_in(new_p) != 0) {
            randPoints.push_back(new_p);
        }
    } while ( std::next_permutation(myints.begin(),myints.end()) );
}*/


// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //

template <class Polytope, class PointList, class Parameters, class Point>
void rand_point_generator(Polytope &P,
                         Point &p,   // a point to start
                         unsigned int rnum,
                         unsigned int walk_len,
                         PointList &randPoints,
                         Parameters &var)  // constants for volume
{
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    unsigned int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    std::vector <NT> lamdas(P.num_of_hyperplanes(), NT(0));
    std::vector <NT> Av(P.num_of_hyperplanes(), NT(0));
    unsigned int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta, lambda;
    Point p_prev = p;
    Point v(n);


    if (var.ball_walk) {
        ball_walk <RNGType> (p, P, ball_rad);
    }else if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else if (var.rdhr_walk) {
        v = get_direction<RNGType, Point, NT>(n);
        std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av);
        var.BoundCalls = var.BoundCalls + 2.0;
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        p = (lambda * v) + p;
    } else {
        billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, true);
    }
    //std::cout<<"is_in = "<<P.is_in(p)<<std::endl;
    //std::cout<<"rnum = "<<rnum<<" walk_step = "<<walk_len<<std::endl;

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk<RNGType>(p, P, ball_rad);
            }else if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas);
            } else if (var.rdhr_walk) {
                v = get_direction<RNGType, Point, NT>(n);
                std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av, lambda);
                var.BoundCalls = var.BoundCalls + 2.0;
                lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
                p = (lambda * v) + p;
            } else {
                billiard_walk(P, p, var.diameter, lamdas, Av, lambda,  var);
            }
        }
        //std::cout<<"is_in = "<<P.is_in(p)<<std::endl;
        //std::cout<<"i = "<<i<<std::endl;
        randPoints.push_back(p);
    }
 
}



template <class BallPoly, class PointList, class Parameters, class Point>
void rand_point_generator(BallPoly &PBLarge,
                         Point &p,   // a point to start
                         unsigned int rnum,
                         unsigned int walk_len,
                         PointList &randPoints,
                         BallPoly &PBSmall,
                         unsigned int &nump_PBSmall,
                         Parameters &var) {  // constants for volume

    typedef typename Point::FT NT;
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    std::vector <NT> lamdas(PBLarge.num_of_hyperplanes(), NT(0));
    std::vector <NT> Av(PBLarge.num_of_hyperplanes(), NT(0));
    unsigned int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta, lambda;
    Point p_prev = p, v(n);

    if (var.ball_walk) {
        ball_walk<RNGType>(p, PBLarge, ball_rad);
    }else if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = PBLarge.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else if (var.rdhr_walk) {
        v = get_direction<RNGType, Point, NT>(n);
        std::pair <NT, NT> bpair = PBLarge.line_intersect(p, v, lamdas, Av);
        var.BoundCalls = var.BoundCalls + 2.0;
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        p = (lambda * v) + p;
    } else {
        billiard_walk(PBLarge, p, var.diameter, lamdas, Av, lambda, var, true);
    }

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk<RNGType>(p, PBLarge, ball_rad);
                //std::cout<<"is in zb_itBW = "<<PBLarge.is_in(p)<<std::endl;
            }else if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, PBLarge, rand_coord, rand_coord_prev, kapa, lamdas);
                //std::cout<<"is in zb_itCD = "<<PBLarge.is_in(p)<<std::endl;
            } else if (var.rdhr_walk) {
                v = get_direction<RNGType, Point, NT>(n);
                std::pair <NT, NT> bpair = PBLarge.line_intersect(p, v, lamdas, Av, lambda);
                var.BoundCalls = var.BoundCalls + 2.0;
                lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
                p = (lambda * v) + p;
            } else {
                billiard_walk(PBLarge, p, var.diameter, lamdas, Av, lambda, var);
            }
        }
        if (PBSmall.second().is_in(p) == -1) {//is in
            randPoints.push_back(p);
            ++nump_PBSmall;
        }
    }
}


template <class Polytope, class Point, class Parameters, typename NT>
void uniform_first_point(Polytope &P,
                               Point &p,   // a point to start
                               Point &p_prev, // previous point
                               unsigned int &coord_prev, // previous coordinate ray
                               unsigned int walk_len, // number of steps for the random walk
                               std::vector<NT> &lamdas,
                               std::vector<NT> &Av,
                               NT &lambda,
                               Parameters &var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n, rand_coord;
    NT kapa;
    Point v(n);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType &rng2 = var.rng;

    if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        coord_prev = rand_coord;
    } else if (var.rdhr_walk) {
        v = get_direction<RNGType, Point, NT>(n);
        std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av);
        var.BoundCalls = var.BoundCalls + 2.0;
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        p = (lambda * v) + p;
    } else {
        billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, true);
    }
    walk_len--;

    for (unsigned int j = 0; j < walk_len; j++) {
        if (var.cdhr_walk) {
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, kapa, lamdas);
            coord_prev = rand_coord;
        } else if (var.rdhr_walk) {
            v = get_direction<RNGType, Point, NT>(n);
            std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av, lambda);
            var.BoundCalls = var.BoundCalls + 2.0;
            lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
            p = (lambda * v) + p;
        } else {
            billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var);
        }
    }
}


template <class Polytope, class Point, class Parameters, typename NT>
void uniform_next_point(Polytope &P,
                        Point &p,   // a point to start
                        Point &p_prev, // previous point
                        unsigned int &coord_prev, // previous coordinate ray
                        unsigned int walk_len, // number of steps for the random walk
                        std::vector<NT> &lamdas,
                        std::vector<NT> &Av,
                        NT &lambda,
                        Parameters &var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n, rand_coord;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType &rng2 = var.rng;
    NT ball_rad = var.delta, kapa;
    Point v(n);


    if (var.ball_walk) {
        for (unsigned int j = 0; j < walk_len; j++) ball_walk<RNGType>(p, P, ball_rad);
    } else if (var.rdhr_walk) {

        for (unsigned int j = 0; j < walk_len; j++) {
            v = get_direction<RNGType, Point, NT>(n);
            std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av, lambda);
            var.BoundCalls = var.BoundCalls + 2.0;
            lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
            p = (lambda * v) + p;
        }
    } else if(var.bill_walk) {
        for (unsigned int j = 0; j < walk_len; j++) billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var);
    }else {
        for (unsigned int j = 0; j < walk_len; j++) {
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, kapa, lamdas);
            coord_prev = rand_coord;
        }
    }
}


// ----- HIT AND RUN FUNCTIONS ------------ //

/*
//hit-and-run with random directions and update
template <class Polytope, class Point, class Parameters>
void hit_and_run(Point &p,
                Polytope &P,
                Parameters &var) {
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, NT>(n);
    std::pair <NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point b1 = (min_plus * l) + p;
    Point b2 = (max_minus * l) + p;
    NT lambda = urdist(rng);
    p = (lambda * b1);
    p = ((1 - lambda) * b2) + p;
}*/


//hit-and-run with orthogonal directions and update
template <class Polytope, class Point, typename NT>
void hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             Polytope &P,
                             unsigned int rand_coord,
                             unsigned int rand_coord_prev,
                             NT kapa,
                             std::vector<NT> &lamdas) {
    std::pair <NT, NT> bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    p_prev = p;
    p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
}


template <class ConvexBody, class Point, class Parameters, typename NT>
void billiard_walk(ConvexBody &P, Point &p, NT che_rad, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev,
        Parameters &var, bool first = false) {

    typedef typename Parameters::RNGType RNGType;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    NT T = urdist(rng) * che_rad;
    Point v = get_direction<RNGType, Point, NT>(n);

    if (first) {

        //std::cout<<"[1]before line intersection"<<std::endl;
        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av);
        var.BoundCalls = var.BoundCalls + 1.0;
        //std::cout<<"[1]after line intersection"<<std::endl;
        if (T <= pbpair.first) {
            p = (T * v) + p;
            lambda_prev = T;
            return;
        }
        lambda_prev = 0.999 * pbpair.first;
        p = (lambda_prev * v) + p;
        T -= lambda_prev;
        //std::cout<<"[1]before reflection"<<std::endl;
        P.compute_reflection(v, p, pbpair.second);
        //std::cout<<"[1]after reflection"<<std::endl;
    }

    while (true) {

        //std::cout<<"[2]before line intersection"<<std::endl;
        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av, lambda_prev);
        var.BoundCalls = var.BoundCalls + 1.0;
        //std::cout<<"[2]after line intersection"<<std::endl;
        //std::cout<<"lambda_pos ="<<pbpair.first<<" T = "<<T<<std::endl;
        if (T <= pbpair.first) {
            p = (T * v) + p;
            lambda_prev = T;
            break;
        }

        lambda_prev = 0.995 * pbpair.first;
        p = (lambda_prev * v) + p;
        //std::cout<<"[in bill] is_in = "<<P.is_in(p)<<std::endl;
        T -= lambda_prev;
        //std::cout<<"[2]before reflection"<<std::endl;
        P.compute_reflection(v, p, pbpair.second);
        //std::cout<<"[2]after reflection"<<std::endl;
    }
}


#endif //RANDOM_SAMPLERS_H
