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

#include "ratio_esti.h"

template <class T>
NT Vpoly_volume (T &P, vars var) {

    int n = var.n, min_index, max_index, index;
    var.delta = 4.0*0.8/NT(n);
    std::vector<Interballs> ConvSet;
    std::vector<Ball> vecBalls;
    P.print();
    bool print = var.verbose;
    NT vol;
    if (print) std::cout<<"\n\ncomputing ball schedule...\n"<<std::endl;
    NT round_value, last_ratio;
    get_ball_schedule(P, ConvSet, vecBalls, 0.70, last_ratio, var);
    if (print) std::cout<<"ball schedule computed!\n"<<std::endl;
    if (print) std::cout<<"number of conv bodies= "<<ConvSet.size()<<std::endl;

    NT rad = vecBalls[0].radius();
    vol = (std::pow(M_PI,n/2.0)*(std::pow(rad, n) ) ) / (tgamma(n/2.0+1));
    std::cout<<"vol of enclosing ball = "<<vol<<std::endl;

    //int mm = ConvSet.size();
    //NT error = 0.1;
    //NT curr_eps = error/std::sqrt((NT(mm)));
    std::vector<NT> ratios;// = ratio_esti(P, ConvSet, 1.0, curr_eps, var);
    std::pair<NT,NT> mv;
    //typename std::vector<NT>::iterator ratIt = ratios.begin();

    //for (; ratIt!=ratios.end(); ++ratIt) {
     //   vol = vol * (*ratIt);
    //}

    if (true) {
    int W = 10*n*n*std::log2(NT(n))+500;
    //int W =100;
    std::vector<NT> last_W2(W,0);
    Point p(n);
    int mm = ConvSet.size();
    std::vector<NT> ratios(mm,0);
    //var.walk_steps=1;
    var.walk_steps = 10+n/10;

    NT error = 0.1, curr_eps, min_val, max_val, val;
    curr_eps = error/std::sqrt((NT(mm)));
    typename  std::vector<Interballs>::iterator CnvIt = ConvSet.begin();
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;

    Interballs S1(n);
    Interballs S2(n);
    int i, count;
    NT countIn, countTot=100.0;
    std::list<Point> randPoints, nextrandPoints;
    nextrandPoints.clear();
    int bef_points;
    for ( i=0;  i<mm-1; i++) {
        S1 = ConvSet[i];
        S2 = ConvSet[i+1];
        p = Point(n);

        countIn = 0.0;
        countTot = 0.0;

        min_val = minNT;
        max_val = maxNT;
        min_index = W-1;
        max_index = W-1;
        index = 0;
        std::vector<NT> last_W=last_W2;
        bool done = false, first_rnum = true;
        count = 0;
        bool reuse = true;
        //nextrandPoints.clear();
        //ratios.clear()

        while(!done){//} && count<60000) {
            randPoints.clear();
            //rand_point(S1, p, var);
            if(reuse){
                randPoints = nextrandPoints;
                bef_points = nextrandPoints.size();
                std::cout<<"reuse = "<<nextrandPoints.size()<<std::endl;
                nextrandPoints.clear();
                reuse = false;
            }else if (first_rnum){
                rand_point_generator(S1, p, ((int)std::pow(1.0,-2.0) * 400 * n * std::log(n))-bef_points, var.walk_steps, randPoints, var);
                std::cout<<"first rnum = "<<((int)std::pow(1.0,-2.0) * 400 * n * std::log(n))<<" done: "<<((int)std::pow(1.0,-2.0) * 400 * n * std::log(n))*((int)10+n/10)<<std::endl;
                first_rnum = false;
            } else {
                done=true;
                break;
                rand_point_generator(S1, p, W, var.walk_steps, randPoints, var);
            }
            //countTot += 1.0;
            count++;
            //countIn = 0.0;

            /*randPoints.clear();
            rand_point_generator(S1, p, 100, var.walk_steps, randPoints, var);
            rpit = randPoints.begin();
            count++;
            for ( ; rpit!=randPoints.end(); ++rpit) {
                if(S2.is_in(*rpit)==-1){
                    countIn += 1.0;
                }
            }*/

            rpit = randPoints.begin();
            for ( ; rpit!=randPoints.end(); ++rpit) {
                countTot += 1.0;

            if (S2.is_in(*rpit)==-1) {
                countIn += 1.0;
                nextrandPoints.push_back(*rpit);
            }
            val = countIn / countTot;
            //ratios.push_back(val);
            //mv = getMeanVariance(ratios);
            //val = mv.first;
            //std::cout<<"val = "<<val<<std::endl;
            //val = countTot / countIn;

            last_W[index] = val;
            if(val<=min_val){
                min_val = val;
                min_index = index;
            }else if(min_index==index){
                minmaxIt = std::min_element(last_W.begin(), last_W.end());
                min_val = *minmaxIt;
                min_index = std::distance(last_W.begin(), minmaxIt);
            }

            if(val>=max_val){
                max_val = val;
                max_index = index;
            }else if(max_index==index){
                minmaxIt = std::max_element(last_W.begin(), last_W.end());
                max_val = *minmaxIt;
                max_index = std::distance(last_W.begin(), minmaxIt);
            }

            //std::cout<<(max_val-min_val)/max_val<<" > "<<curr_eps/2.0<<std::endl;
            if( (max_val-min_val)/max_val<curr_eps/2.0 && false){
                done=true;
                break;
            }

            index = index%W+1;

            if(index==W) index=0;
            }


        }
        if(print) std::cout<<"ratio "<<i<<" = "<<val<<" N_"<<i<<" = "<<countTot<<" num_of_balls = "<<S1.num_of_balls()<<std::endl;
        std::cout<<"walk step = "<<var.walk_steps<<std::endl;
        vol = vol * val;
        //vol = vol * (1.0 / val);

    }

    //var.walk_steps=2*(10+10/n);
    var.walk_steps = 10*std::log2(NT(n));
    //W = W/2;
    S2 = ConvSet[mm-1];
    p = Point(n);

    countTot=0.0;
    countIn = 0.0;
    //countTot = 1200.0;
    //std::cout<<"last_ratio = "<<last_ratio<<std::endl;
    //countIn = last_ratio*1200.0;
    //std::cout<<"countin = "<<countIn<<std::endl;

    W =4*n*n + 500;
    min_val = minNT;
    max_val = maxNT;
    min_index = W-1;
    max_index = W-1;
    index = 0;
    //std::vector<NT> last_W=last_W2;
    std::vector<NT> last_W(W,0);
    bool done = false;
    ratios.clear();


    while(!done) {
        rand_point(S2, p, var);
        countTot += 1.0;

        /*countIn = 0.0;
        randPoints.clear();
        rand_point_generator(S2, p, 100, var.walk_steps, randPoints, var);
        rpit = randPoints.begin();
        for ( ; rpit!=randPoints.end(); ++rpit) {
            if(P.is_in(*rpit)==-1){
                countIn += 1.0;
            }
        }*/


        if (P.is_in(p)==-1) {
            countIn += 1.0;
        }
        val = countIn / countTot;
        //ratios.push_back(val);
        //mv = getMeanVariance(ratios);
        //val = mv.first;
        //val = countTot / countIn;

        last_W[index] = val;
        if(val<=min_val){
            min_val = val;
            min_index = index;
        }else if(min_index==index){
            minmaxIt = std::min_element(last_W.begin(), last_W.end());
            min_val = *minmaxIt;
            min_index = std::distance(last_W.begin(), minmaxIt);
        }

        if(val>=max_val){
            max_val = val;
            max_index = index;
        }else if(max_index==index){
            minmaxIt = std::max_element(last_W.begin(), last_W.end());
            max_val = *minmaxIt;
            max_index = std::distance(last_W.begin(), minmaxIt);
        }

        if( (max_val-min_val)/max_val<=curr_eps / 2.0 ){
            done=true;
        }

        index = index%W+1;

        if(index==W) index=0;


    }
    if(print) std::cout<<"ratio "<<mm-1<<" = "<<val<<" N_"<<mm-1<<" = "<<countTot<<" num_of_balls = "<<S2.num_of_balls()<<std::endl;
    vol = vol * val;
    //vol = vol * (1.0 / val);
    }




    return vol;
}