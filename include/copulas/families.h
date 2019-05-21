// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

#ifndef FAMILIES_H
#define FAMILIES_H


template <class VT, typename NT, class MT>
void get_bodies (MT RetMat, MT allEll, std::vector<MT> &ellipsoids,  MT &hyperplanes, int W) {

    int d = RetMat.cols(), L = RetMat.rows();
    VT hyp(d);// = VT::Ones(d);
    MT compmat(W, d);
    MT temp_ell(d,d);

    for (int i = 0; i < L-W+1; ++i) {
        //compmat = RetMat.block<i+W-1,d-1>(i,0);
        hyp = VT::Ones(d);
        for (int k = 0; k < d; ++k) {
            for (int j = 0; j < W-1; ++j) {
                hyp(k) = hyp(k) * (1.0 + RetMat(i+j, k));
            }
            hyp(k) = hyp(k) - 1.0;
        }
        hyperplanes.row(i) = hyp;

        for (int l = 0; l < d; ++l) {
            temp_ell.row(l) = allEll.row(i*d+l);
        }
        ellipsoids.push_back(temp_ell);
    }
    //allEll.reshape(0,0);
    //RetMat.reshape(0,0);
}



template <class RNGType, class VT, class MT, typename NT>
std::pair<MT, MT> get_constants(std::vector<MT> &ellipsoids,  MT &hyperplanes, std::vector<MT> copulas, int M, int N,
                                std::vector<NT> &mins) {

    MT copula(M, M);
    int d = hyperplanes.row(0).size(), K = ellipsoids.size(), row, col;
    MT points(d, N);
    exp_simplex<NT, RNGType>(d, N, points);

    MT hyp_vals = hyperplanes * points;
    MT hyp_vals2 = hyp_vals;
    MT ell_cons(K, M), hyp_cons(K, M), vecs(d, N), ell_vals(K, N);
    VT row_vals(N), ell_row_vals(N);

    for (int i = 0; i < K; ++i) {
        std::sort(hyp_vals.row(i).data(), hyp_vals.row(i).data() + hyp_vals.row(i).size());
    }

    int count_ell = 0;
    NT pos = 1.0 / NT(M);
    for (typename std::vector<MT>::iterator ellit = ellipsoids.end();  ellit!=ellipsoids.end(); ++ellit, ++count_ell) {
        vecs = (*ellit) * points;
        for (int i = 0; i < N; ++i) {
            ell_vals(count_ell, i) = points.col(i).transpose() * vecs.col(i);
        }
        ell_row_vals = ell_vals.row(count_ell);
        std::sort(ell_vals.row(count_ell).data(), ell_vals.row(count_ell).data()
                                                                            + ell_vals.row(count_ell).size());

        for (int i=1; i<M; i++) {
            row_vals = ell_vals.row(count_ell);
            ell_cons(count_ell, i - 1) = row_vals(((int) std::floor(i * (pos) * (NT(N)))));

            row_vals = hyp_vals.row(count_ell);
            hyp_cons(count_ell, i - 1) = row_vals(((int) std::floor(i * (pos) * (NT(N)))));
        }
        copula = MT::Zero(M,M);
        for (int k = 0; k < N; ++k) {

            row = -1;
            col = 1;
            for (int j = 0; j < M - 1; ++j) {
                if (ell_row_vals(k) < ell_cons(j)) {
                    row = j;
                    break;
                }
            }
            for (int j = 0; j < M - 1; ++j) {
                if (hyp_vals2(count_ell, k) < hyp_cons(j)) {
                    col = j;
                    break;
                }
            }
            if (col == -1) col = M - 1;
            if (row == -1) row = M - 1;
            copula(row, col) = copula(row, col) + 1.0;
        }
        mins.push_back(copula.minCoeff()/NT(N));
        copulas.push_back(copula);
    }

    return std::pair<MT, MT> (ell_cons, hyp_cons);
}

#endif
