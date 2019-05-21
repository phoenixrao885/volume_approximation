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
std::pair<MT, MT> get_constants(std::vector<MT> &ellipsoids,  MT &hyperplanes, std::vector<MT> &copulas, int M, int N,
                                std::vector<NT> &mins) {

    MT copula(M, M);
    int d = hyperplanes.row(0).size(), K = ellipsoids.size(), row, col;
    //std::cout<<"K = "<<K<<" d = "<<d<<std::endl;
    MT points(d, N);
    exp_simplex<NT, RNGType>(d, N, points);
    //std::cout<<"points are sampled!"<<std::endl;
    std::cout<<points.cols()<<" "<<points.rows()<<"\n"<<std::endl;


    MT hyp_vals = hyperplanes * points;
    std::cout<<hyp_vals.cols()<<" "<<hyp_vals.rows()<<"\n"<<std::endl;
    MT hyp_vals2 = hyp_vals;
    MT ell_cons(K, M-1), hyp_cons(K, M-1), vecs(d, N), ell_vals(K, N);
    VT row_vals(N), ell_row_vals(N);

    for (int i = 0; i < K; ++i) {
        row_vals = hyp_vals.row(i);
        std::sort(row_vals.data(), row_vals.data() + row_vals.size());
        hyp_vals.row(i) = row_vals;
        //std::cout<<hyp_vals.row(i)<<"\n"<<std::endl;
    }
    //std::cout<<"hello0.1"<<std::endl;
    //std::cout<<hyp_vals<<"\n"<<std::endl;

    int count_ell = 0;
    NT pos = 1.0 / NT(M);
    for (typename std::vector<MT>::iterator ellit = ellipsoids.begin();  ellit!=ellipsoids.end(); ++ellit, ++count_ell) {
        //std::cout<<"hello0.3"<<std::endl;
        vecs = (*ellit) * points;
        //std::cout<<"hello0.5"<<std::endl;
        for (int i = 0; i < N; ++i) {
            //std::cout<<"hello1"<<std::endl;
            //ell_vals(count_ell, i) = points.col(i).transpose() * vecs.col(i);
            ell_vals(count_ell, i) = points.col(i).dot(vecs.col(i));
            //std::cout<<"hello2"<<std::endl;
        }
        ell_row_vals = ell_vals.row(count_ell);
        row_vals = ell_vals.row(count_ell);
        std::sort(row_vals.data(), row_vals.data() + row_vals.size());
        ell_vals.row(count_ell) = row_vals;

        row_vals = ell_vals.row(count_ell);
        for (int i=1; i<M; i++) {
            ell_cons(count_ell, i - 1) = row_vals(((int) std::floor(i * (pos) * (NT(N)))));
        }

        row_vals = hyp_vals.row(count_ell);
        for (int i=1; i<M; i++) {
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
        std::cout<<"number of copulas = "<<copulas.size()<<std::endl;
    }

    return std::pair<MT, MT> (ell_cons, hyp_cons);
}

#endif
