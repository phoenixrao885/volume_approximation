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
NT Vpoly_volume (T &P, vars var) {

    std::vector<InterBalls> ConvSet;
    int n = var.n;
    bool print = var.verbose;
    NT vol;

    if (print) std::cout<<"Computing minimum enclosing ball..\n";
    std::pair<Point, NT> res = compute_minball(P);
    if (print) std::cout<<"Computation of minimum enclosing ball completed!\n";
    Point xc = res.first;
    if (print) std::cout<<"center of minimum enclosing ball = ";
    xc.print();
    if (print) std::cout<<"radius of minimum enclosing ball = "<<res.second<<std::endl;

    if (print) std::cout<<"\nshifting polytope...\n\n";
    VT c_e(n);
    for(int i=0; i<n; i++){
        c_e(i)=xc[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);

    xc = Point(n);

    std::vector<ball> S0;
    S0.push_back(ball(c, res.second * res.second));
    // starting convex body is the minimum enclosing ball of the V-polytope
    ConvSet.push_back(P.dimension(), InterBalls(S0));

    vol = (std::pow(M_PI,n/2.0)*(res.second, n) )  / (tgamma(n/2.0+1));

    return vol;
}