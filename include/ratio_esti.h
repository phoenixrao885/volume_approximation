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


template <class T1, class T2, typename FT>
std::vector<FT> ratio_esti(T1 &P, std::vector<T2> ConvSet, FT a, FT err, vars &var){

    int mm = ConvSet.size(), n =var.n, ni = 0;
    bool print = var.verbose;
    var.walk_steps = 1;
    T2 S1(n);
    T2 S2(n);
    Point p(n);
    std::list<Point> randPoints;
    typename std::list<Point>::iterator rpit;
    bool done;
    FT totcount = 300, countIn = 0.0, val, meanval, dx, stdev, meanvalbef;
    std::vector<FT> ratios, retratios;
    std::pair<FT,FT> mv;
    int contpos = 0;

    for (int i = 0; i < mm-1; ++i) {
        if (i==mm-2){
            var.walk_steps=2*(10+10/n);
        }
        S1 = ConvSet[i];
        S2 = ConvSet[i+1];
        p = Point(n);
        done = false;
        ni = 0;
        ratios.clear();
        meanvalbef = 0.0;

        while (!done) {
            randPoints.clear();

            if (i==mm-2) {
                rand_point_generator(S2, p, 300, var.walk_steps, randPoints, var);
            } else {
                rand_point_generator(S1, p, 300, var.walk_steps, randPoints, var);
            }
            countIn = 0.0;

            rpit = randPoints.begin();
            if (i==mm-2) {
                for ( ; rpit!=randPoints.end(); ++rpit) {
                    if(P.is_in(*rpit)==-1){
                        countIn += 1.0;
                    }
                }
            } else {
                for ( ; rpit!=randPoints.end(); ++rpit) {
                    if(S2.is_in(*rpit)==-1){
                        countIn += 1.0;
                    }
                }
            }

            ratios.push_back(countIn / totcount);
            ni++;

            if(ni>=10) {
                mv = getMeanVariance(ratios);
                meanval = mv.first;
                //std::cout<<"meanval = "<<meanval<<" val_bef = "<<meanvalbef<<" err = "<<std::abs(meanval - meanvalbef) / meanval<<" dem err = "<<err/2.0<<" ni = "<<ni<<std::endl;
                if ((std::abs(meanval - meanvalbef)) / meanval < err / 2.0) {
                    contpos++;
                    if (contpos >= 20) {
                        retratios.push_back(meanval);
                        if(print) std::cout<<"ratio "<<i<<" = "<<meanval<<" N_"<<i<<" = "<<300*ni<<std::endl;
                        done = true;
                    }
                    //contpos++;
                } else {
                    contpos = 0;
                }
                meanvalbef = meanval;
            }
        }
        /*
                stdev = std::sqrt(mv.second);
                dx = stdev;
                std::cout<<"err = "<<err<<" var = "<<mv.second<<" std = "<<stdev<<" dx = "<<dx<<" ni = "<<ni<<std::endl;
                std::cout<<"meanval = "<<meanval<<" t-value = "<<((1 + err)/err)*dx + (0.7*FT(ni-1)*stdev)/std::sqrt(FT(ni))<<" ni = "<<ni<<std::endl;
                if ( meanval > ((1 + err)/err)*dx + (0.7*FT(ni-1)*stdev)/std::sqrt(FT(ni)) ) {
                    std::cout<<"meanval = "<<meanval<<" ni = "<<ni<<std::endl;
                    contpos++;
                    if (contpos==10){
                        retratios.push_back(meanval);
                        if(print) std::cout<<"ratio "<<i<<" = "<<meanval<<" N_"<<i<<" = "<<120*ni<<std::endl;
                        done = true;
                    }
                } else {
                    contpos=0;
                }
            }

        }*/



    }
    return retratios;

}