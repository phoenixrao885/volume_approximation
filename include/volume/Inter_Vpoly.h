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

#ifndef INTER_VPOLY_H
#define INTER_VPOLY_H

template <class Ball, class HPolytope, class VPolytope, typename NT, class Parameters>
NT Inter_Vpoly_Vol(VPolytope &P, NT &p_value, NT epsilon, bool round, Parameters &var) {

    typedef typename VPolytope::PolytopePoint Point;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;
    typedef BallIntersectPolytope<HPolytope, Ball> HPolyBall;
    //typedef Ball<Point> Ball;
    //typedef IntersectionOfBalls<Ball, NT> Interballs;
    NT maxNT = 1.79769e+308;
    NT minNT = -1.79769e+308;
    int n = var.n, min_index, max_index, index;
    //var.delta = 4.0*0.8/NT(n);
    std::vector<HPolyBall> ConvSet;
    std::vector<Ball> vecBalls;
    P.print();
    bool print = var.verbose;
    NT vol;
    NT round_value = 1.0;

    /*
    if (round) {
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
        VT c_e(n);
        for (int l = 0; l < n; ++l) {
            c_e(l) = xc[l];
        }
        P.shift(c_e);
        std::pair<Point,NT> Che(Point(n), 0.1);
        std::pair <NT, NT> res_round;
        res_round = rounding_min_ellipsoid(P, Che, var);
        round_value = res_round.first;
    }*/



    if (print) std::cout<<"\n\ncomputing schedule...\n"<<std::endl;
    NT last_ratio, aaa= 0.7;
    bool empty = false;
    Point xc = P.getInnerPoint(empty);
    if (empty) {
        std::cout<<"intersection is empty"<<std::endl;
        return -1.0;
    }
    get_hyperplane_annealing2(P, ConvSet, p_value, aaa, epsilon, xc, var);
    //if (print) std::cout<<"ball schedule computed!\n"<<std::endl;
    if (print) std::cout<<"number of conv bodies= "<<ConvSet.size()<<std::endl;
    if (print) std::cout<<"number of hyperplanes in last convex = "<<ConvSet[ConvSet.size()-1].num_of_hyperplanes()<<std::endl;
    NT rad = ConvSet[0].second().radius();
    vol = (std::pow(M_PI,n/2.0)*(std::pow(rad, NT(n)) ) ) / (tgamma(n/2.0+1));
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
        var.coordinate=true;

        NT error = 0.1, curr_eps;
        curr_eps = error/std::sqrt((NT(mm)));
        typename  std::vector<HPolyBall>::iterator CnvIt = ConvSet.begin();
        typename std::vector<NT>::iterator minmaxIt;
        typename std::list<Point>::iterator rpit;

        HPolyBall S1;
        HPolyBall S2;
        int i, count;
        NT countIn, countTot=100.0, min_val, max_val, val;
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
                    rand_point_generator(S1, p, (2*(int)std::pow(1.0,-2.0) * 400 * n * std::log(n))-bef_points, var.walk_steps, randPoints, var);
                    std::cout<<"first rnum = "<<(2*(int)std::pow(1.0,-2.0) * 400 * n * std::log(n))<<" done: "<<((int)std::pow(1.0,-2.0) * 400 * n * std::log(n))*((int)10+n/10)<<std::endl;
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
            if(print) std::cout<<"ratio "<<i<<" = "<<val<<" N_"<<i<<" = "<<countTot<<" num_of_hyperplanes = "<<S1.num_of_hyperplanes()<<std::endl;
            std::cout<<"walk step = "<<var.walk_steps<<std::endl;
            vol = vol * val;
            //vol = vol * (1.0 / val);

        }

        //var.walk_steps=2*(10+10/n);
        var.walk_steps = 20*std::log2(NT(n));
        //W = W/2;
        S2 = ConvSet[mm-1];
        p = Point(n);

        countTot=0.0;
        countIn = 0.0;
        //countTot = 1200.0;
        //std::cout<<"last_ratio = "<<last_ratio<<std::endl;
        //countIn = last_ratio*1200.0;
        //std::cout<<"countin = "<<countIn<<std::endl;

        std::cout<<"dimension = "<<n<<"W = "<<var.walk_steps<<std::endl;
        W =10*n*n + 1200;
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
        if(print) std::cout<<"ratio "<<mm-1<<" = "<<val<<" N_"<<mm-1<<" = "<<countTot<<" num_of_hyperplanes = "<<S2.num_of_hyperplanes()<<std::endl;
        vol = vol * val;
        //vol = vol * (1.0 / val);
    }




    return round_value * vol;
}

#endif
