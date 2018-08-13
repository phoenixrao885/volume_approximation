// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//
// For the computation of the minimum volume enclosing ball of the V-polytope we use the Seb library
//
// Authors: Martin Kutz <kutz@math.fu-berlin.de>,
//          Kaspar Fischer <kf@iaeth.ch>

#include <iostream>
#include <cstdio>

#include "Seb.h"


template <class T>
std::pair<Point, NT> compute_minball(T &P) {

    typedef Seb::Point<NT> Point_seb;
    typedef std::vector<Point_seb> PointVector;
    typedef Seb::Smallest_enclosing_ball<NT> Miniball;

    int d = P.dimension(), i ,j;
    int k = P.num_of_vertices();

    PointVector S;
    std::vector<NT> coords(d);
    MT V = P.get_mat();

    for (i=0; i<k; i++) {
        for (j=0; j<d; j++) {
            coords[j] = V(i, j);
        }
        S.push_back(Point_seb(d,coords.begin()));
    }

    // Compute the miniball by inserting each value
    Miniball mb(d, S);

    NT rad = mb.radius();

    Miniball::Coordinate_iterator center_it = mb.center_begin();
    for (j=0; j<d; ++j)
        coords[j] =  center_it[j];

    Point xc(d, coords.begin(), coords.end());

    return std::pair<Point, NT> (xc, rad);

}
