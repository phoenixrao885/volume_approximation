// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef EST_RATIO_ONE_H
#define EST_RATIO_ONE_H

template <class Zonotope, class HPolytope, typename NT, class Parameters>
NT est_ratio_hzono(Zonotope &Z, HPolytope &HP, NT error, Parameters &var, NT &steps) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    int W=4*n*n+500;
    int m = Z.num_of_generators();
    NT curr_eps = error;
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int min_steps=0;
    std::vector<NT> last_W(W,0);
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    NT val;

    // Set the radius for the ball walk if it is requested
    //if (var.ball_walk) {
    //if (var.deltaset) {
    //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
    //}
    //}

    //MT sigma2;
    //MT sample;
    NT countIn = 0.0;
    NT totCount = 0.0;
    Point p(n);
    while(!done){

        //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
        //sigma2 = (1.0/(2.0*(*avalsIt)))*sigma;
        rand_point(HP, p, var);
        if(Z.is_in(p)==-1) {
            countIn = countIn + 1.0;
            //q2=*rpit;
        }
        totCount = totCount + 1.0;

        //pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);

        //for (int k = 0; k < 4*W; ++k) {
        //*itsIt = *itsIt + 1.0;
        //*fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
        // *fnIt = *fnIt + std::exp(-(*(avalsIt + 1))*(pointset.col(k).squaredNorm())) / std::exp(-(*avalsIt)*(pointset.col(k).squaredNorm()));
        //val = (*fnIt) / (*itsIt);

        val = countIn / totCount;
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

        if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
            std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
            steps = (totCount - 1200.0);
            return val;
        }

        index = index%W+1;

        if(index==W) index=0;

    }
    std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
    return val;
}

template <class Zonotope, class MT, class VT, typename NT>
NT est_ratio_zono(Zonotope &Z, NT &prob, NT error, int W, Rcpp::Function rtmvnorm,
                    Rcpp::Function mvrandn, MT &sigma, MT &G, VT &l, VT &u) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    W=4*n*n+500;
    int m = Z.num_of_generators();
    NT curr_eps = error;
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int min_steps=0;
    std::vector<NT> last_W(W,0);
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    NT val;

    // Set the radius for the ball walk if it is requested
    //if (var.ball_walk) {
    //if (var.deltaset) {
    //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
    //}
    //}

    MT sigma2;
    MT sample;
    NT countIn = 0.0;
    NT totCount = 0.0;
    while(!done){

        //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
        //sigma2 = (1.0/(2.0*(*avalsIt)))*sigma;
        if(prob>0.001) {
            sample = sampleTr(l, u , sigma, 3*W, mvrandn, G);
        } else {
            sample = sampleTr_gibbs(l, u, sigma, 3*W, m*m/10, rtmvnorm, G);
        }
        randPoints.clear();
        for (int i = 0; i < 3*W; ++i) {
            Point p(n, typename std::vector<NT>::iterator(sample.col(i).data()), typename std::vector<NT>::iterator(sample.col(i).data() + n));
            randPoints.push_back(p);
        }
        rpit = randPoints.begin();

        //Point q2;
        for ( ;  rpit!=randPoints.end(); ++rpit) {
            if(Z.is_in(*rpit)==-1) {
                countIn = countIn + 1.0;
                //q2=*rpit;
            }
            totCount = totCount + 1.0;

        //pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);

        //for (int k = 0; k < 4*W; ++k) {
            //*itsIt = *itsIt + 1.0;
            //*fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
           // *fnIt = *fnIt + std::exp(-(*(avalsIt + 1))*(pointset.col(k).squaredNorm())) / std::exp(-(*avalsIt)*(pointset.col(k).squaredNorm()));
            //val = (*fnIt) / (*itsIt);

            val = countIn / totCount;
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

            if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
                std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
                done=true;
                return val;
            }

            index = index%W+1;

            if(index==W) index=0;
        }
    }
    std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
    return val;
}


template <class Point, class Zonoball, class MT, class VT, typename NT, class Parameters>
NT est_ratio_zball_sym(Zonoball ZB, MT sigma, MT G, MT Q0, VT l, VT u, NT delta, NT ratio, NT error, Parameters &var){

    NT countIn = ratio*NT(1200);
    NT totCount = NT(1200);

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    //typedef typename Zonotope::PolytopePoint Point;
    int n = var.n;
    int W=4*n*n+500;
    int m = G.cols();
    NT curr_eps = error;
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int min_steps=0;
    std::vector<NT> last_W(W,0);
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    NT val = ratio;

    // Set the radius for the ball walk if it is requested
    //if (var.ball_walk) {
    //if (var.deltaset) {
    //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
    //}
    //}

    MT sigma2;
    MT sample;
   // NT countIn = 0.0;
   // NT totCount = 0.0;
    Point p(n);
    //std::cout<<"CountIn = "<<countIn<<" totCount = "<<totCount<<std::endl;
    while(!done){

        rand_point(ZB, p, var);

        //Point q2;
        //for ( ;  rpit!=randPoints.end(); ++rpit) {
        if (is_in_sym2(p, Q0, G, delta)) {
            countIn = countIn + 1.0;
        }
            totCount = totCount + 1.0;
        //std::cout<<"CountIn = "<<countIn<<" totCount = "<<totCount<<std::endl;

            //pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);

            //for (int k = 0; k < 4*W; ++k) {
            //*itsIt = *itsIt + 1.0;
            //*fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
            // *fnIt = *fnIt + std::exp(-(*(avalsIt + 1))*(pointset.col(k).squaredNorm())) / std::exp(-(*avalsIt)*(pointset.col(k).squaredNorm()));
            //val = (*fnIt) / (*itsIt);

            val = countIn / totCount;
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

            if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
                std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
                done=true;
                return val;
            }

            index = index%W+1;

            if(index==W) index=0;
        //}
    }

}


template <class Point, class convexB, class Zonotope, class Parameters, typename NT>
NT est_ratio_zonoballs(Zonotope &Z, convexB &b1, NT ratio, NT error, Parameters &var, NT &steps){

    NT countIn = ratio*NT(1200);
    NT totCount = NT(1200);
    //typedef typename Zonotope::PolytopePoint Point;
    //ball b1 = zb1.second();
   // ball b2 = zb2.second();

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    //typedef typename Zonotope::PolytopePoint Point;
    int n = var.n;
    int W=4*n*n+500;
    //int m = Z.num_of_generators();
    NT curr_eps = error;
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int min_steps=0;
    std::vector<NT> last_W(W,0);
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    NT val = ratio;

    // Set the radius for the ball walk if it is requested
    //if (var.ball_walk) {
    //if (var.deltaset) {
    //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
    //}
    //}

    //NT countIn = 0.0;
    //NT totCount = 0.0;
    Point p(n);
    while(!done){

        rand_point(Z, p, var);
        //Point q2;
        //for ( ;  rpit!=randPoints.end(); ++rpit) {
        if (b1.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }

        totCount = totCount + 1.0;

        //pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);

        //for (int k = 0; k < 4*W; ++k) {
        //*itsIt = *itsIt + 1.0;
        //*fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
        // *fnIt = *fnIt + std::exp(-(*(avalsIt + 1))*(pointset.col(k).squaredNorm())) / std::exp(-(*avalsIt)*(pointset.col(k).squaredNorm()));
        //val = (*fnIt) / (*itsIt);

        val = countIn / totCount;
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

        if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
            std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
            steps = (totCount - 1200.0);
            return val;
        }

        index = index%W+1;

        if(index==W) index=0;
        //}
    }

}


#endif