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

template <class T1, typename FT>
void get_next_schedule(T1 &P, std::vector<Interballs> &ConvSet, FT a, vars &var, bool &done){
    Interballs Si = (*ConvSet.end());

}


template <class T1, typename FT>
void get_ball_schedule(T1 &P, std::vector<Interballs> &ConvSet, FT a, vars &var) {

    bool done = false;
    while(!done) {
        get_next_convex(P, ConvSet, a, var, done);
    }
}