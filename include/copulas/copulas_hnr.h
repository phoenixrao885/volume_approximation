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
    int tot_points, count_ell, row, col, count = N;

    tot_points = int( ((1.0+error)/error)*((1.0+error)/error)*zp*zp*((1.0-min_ratio)/min_ratio) );

    typename std::vector<MT>::iterator ellit;
    MT points(d, N), cons(K,N), vecs(d, N);

    while (count < tot_points) {

        exp_simplex<NT, RNGType>(d, N, points);
        count += N;
        ellit = ellipsoids.begin();

        cons = hyperplanes * points;
        count_ell = 0;

        for ( ;  ellit!=ellipsoids.end(); ++ellit, ++count_ell) {

            copula = copulas[count_ell];
            vecs = (*ellit) * points;
            hyp_vals = cons.row(count_ell);
            for (int i = 0; i < N-1; ++i) {

                val_ell = points.col(i).transpose() * vecs.col(i);
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


#endif
