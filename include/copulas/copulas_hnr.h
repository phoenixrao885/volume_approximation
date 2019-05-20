// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

#ifndef COPULAS_HNR_H
#define COPULAS_HNR_H


template <class MT, class VT, typename NT>
std::vector<MT> get_copulas_uniform (std::vector<MT> &ellipsoids, MT &hyperplanes, NT error, NT prob, int M, int N) {

    int d = hyperplanes.row(0).size(), K = ellipsoids.size();
    NT min, val_ell, val_hyp;
    std::vector<NT> mins;
    VT ell_consts(M), hyp_cosnts(M), hyp_vals(N);
    MT copula(M,M);
    std::vector<MT> copulas;
    std::pair<MT, MT> constants = get_constants(ellipsoids, hyperplanes, copulas, M, N, mins);
    ells_consts = constants.first;
    hyps_consts = constants.second;

    int tot_points, count_ell, row, col;
    count = N;
    typename std::vector<MT>::iterator ellit;
    MT points(d, N), cons(K,N), ells_consts(K,M), hyps_consts(K,M), vecs(d, N);

    while (count < tot_points) {

        exp_simplex(d, N, points);
        count += N;
        ellit = ellipsoids.begin();

        cons = hyperplanes * points;
        count_ell = 0;

        for ( ;  ellit!=ellipsoids.end(); ++ellit, ++count_ell) {

            copula = copulas(count_ell)
            vecs = (*ellit) * points;
            hyp_vals = cons.row(count_ell)
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
                    if (val_hyp < hyp_cosnts(j)){
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
