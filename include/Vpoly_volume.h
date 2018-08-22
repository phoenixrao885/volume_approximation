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

    int n = var.n, min_index, max_index, index;
    var.delta = 4.0*0.8/NT(n);
    std::vector<Interballs> ConvSet;
    std::vector<Ball> vecBalls;
    P.print();
    bool print = var.verbose;
    NT vol;
    if (print) std::cout<<"\n\ncomputing ball schedule...\n"<<std::endl;
    NT round_value;
    get_ball_schedule(P, ConvSet, vecBalls, 0.70, var);
    if (print) std::cout<<"ball schedule computed!\n"<<std::endl;
    if (print) std::cout<<"number of conv bodies= "<<ConvSet.size()<<std::endl;

    NT rad = vecBalls[0].radius();
    vol = (std::pow(M_PI,n/2.0)*(std::pow(rad, n) ) ) / (tgamma(n/2.0+1));

    int W = 2*n*n*((int)std::log2(n))+800;
    //int W =1200;
    std::vector<NT> last_W2(W,0);
    Point p(n);
    int mm = ConvSet.size();
    std::vector<NT> ratios(mm,0);
    //var.walk_steps=1;
    var.walk_steps = 1;

    NT error = 0.1, curr_eps, min_val, max_val, val;
    curr_eps = error/std::sqrt((NT(mm)));
    typename  std::vector<Interballs>::iterator CnvIt = ConvSet.begin();
    typename std::vector<NT>::iterator minmaxIt;

    Interballs S1(n);
    Interballs S2(n);
    int i;
    NT countIn, countTot;
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
        bool done = false;


        while(!done) {
            rand_point(S1, p, var);
            countTot += 1.0;

            if (S2.is_in(p)==-1) {
                countIn += 1.0;
            }
            val = countIn / countTot;
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

            if( (max_val-min_val)/max_val<curr_eps/2.0){
                done=true;
            }

            index = index%W+1;

            if(index==W) index=0;


        }
        if(print) std::cout<<"ratio "<<i<<" = "<<val<<" N_"<<i<<" = "<<countTot<<std::endl;
        vol = vol * val;
        //vol = vol * (1.0 / val);

    }

    var.walk_steps=2*(10+10/n);
    //W = W/2;
    S2 = ConvSet[mm-1];
    p = Point(n);

    countIn = 0.0;
    countTot = 0.0;

    min_val = minNT;
    max_val = maxNT;
    min_index = W-1;
    max_index = W-1;
    index = 0;
    //std::vector<NT> last_W=last_W2;
    std::vector<NT> last_W(W,0);
    bool done = false;


    while(!done) {
        rand_point(S2, p, var);
        countTot += 1.0;

        if (P.is_in(p)==-1) {
            countIn += 1.0;
        }
        val = countIn / countTot;
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
    if(print) std::cout<<"ratio "<<mm-1<<" = "<<val<<" N_"<<mm-1<<" = "<<countTot<<std::endl;
    vol = vol * val;
    //vol = vol * (1.0 / val);




    return vol;
}