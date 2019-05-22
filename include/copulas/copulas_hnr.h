// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

#ifndef COPULAS_HNR_H
#define COPULAS_HNR_H


template < class RNGType, class VT, class MT, typename NT>
std::vector<MT> get_copulas_uniform (std::vector<MT> &ellipsoids, MT &hyperplanes, NT error, NT prob, int M, int N) {

    int d = hyperplanes.row(0).size(), K = ellipsoids.size();
    NT min, val_ell, val_hyp;
    std::vector<NT> mins;
    VT ell_consts(M), hyp_consts(M), hyp_vals(N);
    MT copula(M,M);
    std::vector<MT> copulas;
    std::pair<MT, MT> constants = get_constants<RNGType, VT>(ellipsoids, hyperplanes, copulas, M, N, mins);
    MT ells_consts = constants.first;
    MT hyps_consts = constants.second;

    typename std::vector<NT>::iterator minIt;
    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0));
    minIt = std::min_element(mins.begin(), mins.end());
    NT min_ratio = *minIt;
    if (min_ratio==0.0) min_ratio = (1.0/NT(M))*(1.0/NT(M));
    int tot_points, count_ell, row, col, count = N;

    std::cout<<"error = "<<error<<" prob = "<<prob<<" zp = "<<zp<<" min_ratio = "<<min_ratio<<std::endl;
    tot_points = int( ((1.0+error)/error)*((1.0+error)/error)*zp*zp*((1.0-min_ratio)/min_ratio) );
    tot_points=1000000;
    std::cout<<"totpoints = "<<tot_points<<std::endl;
    std::cout<<"number of copulas = "<<copulas.size()<<std::endl;

    typename std::vector<MT>::iterator ellit;
    MT points(d, N), cons(K,N), vecs(d, N);

    while (count < tot_points) {
        std::cout<<"count = "<<count<<std::endl;

        exp_simplex<NT, RNGType>(d, N, points);
        count += N;
        ellit = ellipsoids.begin();

        cons = hyperplanes * points;
        count_ell = 0;

        for ( ;  ellit!=ellipsoids.end(); ++ellit, ++count_ell) {

            copula = MT::Zero(M,M);
            vecs = (*ellit) * points;
            hyp_vals = cons.row(count_ell);
            //std::cout<<"hello"<<std::endl;
            for (int i = 0; i < N; ++i) {

                val_ell = points.col(i).dot(vecs.col(i));
                val_hyp = hyp_vals[i];

                ell_consts = ells_consts.row(count_ell);
                hyp_consts = hyps_consts.row(count_ell);

                row = -1;
                col = -1;
                for (int j = 0; j < M-1; ++j) {
                    if (val_ell < ell_consts(j)){
                        row = j;
                        break;
                    }
                }
                for (int j = 0; j < M-1; ++j) {
                    if (val_hyp < hyp_consts(j)){
                        col = j;
                        break;
                    }
                }
                if (col == -1) col = M - 1;
                if (row == -1) row = M - 1;
                copula(row, col) = copula(row, col) + 1.0;

            }
            copulas[count_ell] = copulas[count_ell] + copula;
        }


    }
    for (int k = 0; k < K; ++k) {
        copulas[k] = copulas[k] * 1.0/NT(count);
    }

    return copulas;

}

/*
template < class RNGType, class VT, class MT, typename NT>
std::vector<MT> get_copulas_hnr (std::vector<MT> &ellipsoids, MT &hyperplanes, NT error, NT prob, int M, int N) {

    int d = hyperplanes.row(0).size(), K = ellipsoids.size();
    NT min, val_ell, val_hyp, lambda, lambda_prev;
    std::vector<NT> mins, hyp_ratios;
    VT ell_consts(M), hyp_consts(M), hyp_vals(N);
    std::vector<MT> copulas;
    VT p(d), p_prev(d);
    std::vector<int> crow, ccol;
    std::pair<MT, MT> constants = get_constants<RNGType, VT>(ellipsoids, hyperplanes, copulas, M, N, mins);
    MT ells_consts = constants.first, hyps_consts = constants.second, copula(M,M), T = Zero(d,d);

    for (int i = 0; i < d; ++i) {
        if(i==d-1) {
            for (int j = 0; j < d; ++j) {
                T(i, j) = -1.0;
            }
        } else {
            T(i,i) = 1.0;
        }
    }
    p(d-1) = 0.0;
    MT temp_mat(d,d), D(K,d), V(K,d), Js(K,d);
    VT Vals(K), dVals(K), Cvals(K), dCvals(K);
    for (int k = 0; k < K; ++k) {
        ells_consts.row(k) = ells_consts.row(k) - VT::Ones(M-1)*ellipsoids[k](d-1,d-1);
        hyps_consts.row(k) = hyps_consts.row(k) - VT::Ones(M-1)*hyperplanes(k,d-1);
        hyperplanes = hyperplanes * T;
        temp_mat = ellipsoids[k]*T;
        Js.row(k) = 2.0 * temp_mat.row(d-1);
        temp_mat = T.transpose()*ellipsoids[k]*T;
        D.row(k) = temp_mat.diagonal();
        V.row(k) = temp_mat*p;
        Vals(k) = p.dot(V.row(k));
        ellipsoids = T.transpose()*ellipsoids[k]*T;
    }

    typename std::vector<MT>::iterator ellit;
    MT points(d, N), cons(K,N), vecs(d, N);

    while (count < tot_points) {

        coord;
        lambda = update_coord();
        lambda_prev = lambda;
        dVals = 2.0*lambda*V.col(coord) + lambda*lambda*D.col(coord) + lambda * Js.col(coord);
        Vals = Vals + dVals;
        dCvals = lambda*hyperplanes.col(coord);

        for (int i = 0; i < K; ++i) {

            col = ccol(i);
            row = crow(i);
            if (dVals(i) >0.0) {
                for (int j = row; j < M - 1; ++j) {
                    if (Vals(i) < ells_consts(i, j)) {
                        crow(i) = j;
                        break;
                    }
                }
                if (ccol(i) == col) ccol(i) = M - 1;
            } else {
                for (int j = 0; j <=row; ++j) {
                    if (Vals(i) < ells_consts(i, j)) {
                        crow(i) = j;
                        break;
                    }
                }
            }
            if (dCvals(i)>0) {
                for (int j = col; j < M-1; ++j) {
                    if (Cvals(i) < hyps_consts(i, j)){
                        ccol(i) = j;
                        break;
                    }
                }
                if (ccol(i) == col) ccol(i) = M - 1;

            } else {
                for (int j = 0; j <=col; ++j) {
                    if (Cvals(i) < hyps_consts(i, j)){
                        ccol(i) = j;
                        break;
                    }
                }
            }
            //copula(crow(i), ccol(i)) = copula(crow(i), ccol(i)) + 1.0;
            copulas[i] = copulas[i] + copula;
            V.row(i) = V.row(i) + ellipsoids[i].col(coord);
        }
    }

}*/


#endif
