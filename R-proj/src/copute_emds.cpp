

#include <Rcpp.h>
#include <RcppEigen.h>
#undef Realloc
#undef Free
#include "lp_lib.h"


template <class MT, class VT, typename NT>
NT get_dis(MT A, VT b, int m, VT f, NT *row, int *colno) {

    lprec *lp;
    lp = make_lp(2 * m + 1, m * m);

    set_add_rowmode(lp, TRUE);

    REAL infinite = get_infinite(lp);
    for (int j = 0; j < m * m; j++) {
        colno[j] = j + 1;
        row[j] = f(j);
        set_bounds(lp, j + 1, 0.0, infinite);
    }
    set_obj_fnex(lp, m * m, row, colno);


    for (int j = 0; j < m * m; j++) {
        row[j] = 1.0; //Aeq(i,j);
    }
    add_constraintex(lp, m * m, row, colno, EQ, 1.0); //beq(i)

    for (int i = 0; i < 2 * m; ++i) {
        for (int j = 0; j < m * m; j++) {
            row[j] = A(i, j); //Aeq(i,j);
        }

        add_constraintex(lp, m * m, row, colno, LE, b(i));
    }

    set_add_rowmode(lp, FALSE);

    set_minim(lp);
    set_verbose(lp, NEUTRAL);

    if (solve(lp) != OPTIMAL) std::cout << "OPA" << std::endl;

    NT res = NT(get_objective(lp));
    delete_lp(lp);

    return res;

}


// [[Rcpp::export]]
Rcpp::NumericVector emd_Bmat (Rcpp::NumericMatrix B, Rcpp::NumericMatrix A, int m,  Rcpp::NumericVector fin, int i=0) {

    typedef double NT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    MT BB = Rcpp::as<MT>(B);
    //std::cout<<B<<std::endl;

    //MT A1(m, m * m);// = matrix(0, m, m * n)
    //MT A2(m, m * m);// = matrix(0, n, m * n)

    int k;
    //for (int i = 0; i < m; i++) {
    //    for (int j = 0; j < m; j++) {
    //        k = j + i * m;
    //        A1(i, k) = 1;
    //        A2(j, k) = 1;
    //    }
    //}

    //MT A(2 * m, m * m);
    //for (int i = 0; i < m; ++i) A.row(i) = A1.row(i);
    //for (int i = m; i < 2*m; ++i) A.row(i) = A2.row(i-m);

    MT AA = Rcpp::as<MT>(A);

    //MT Aeq = MT::Ones(2 * m, m * m);
    //VT beq = VT::Ones(2 * m);
    //VT l = VT::Zero(m * m);
    VT f = Rcpp::as<VT>(fin);
    //std::cout<<AA.rows()<<A.cols()<<std::endl;


    int *colno = NULL;
    REAL *row = NULL;

    colno = (int *) malloc((m * m) * sizeof(*colno));
    row = (REAL *) malloc((m * m) * sizeof(*row));

    int K = BB.cols();
    Rcpp::NumericVector res(K);
    VT bb(2*m);
    for (int k = 0; k < K; ++k) {

      bb = BB.col(k);
      std::cout<<"i = "<<i<<" j = "<<k+i+1<<std::endl;
      res(k) = NT(get_dis(AA, bb, m, f, row, colno));

    }

    free(row);
    free(colno);

    return res;
}

/*
    lprec *lp, *lpcopy;
    lp = make_lp(2 * m + 1, m * m);

    colno = (int *) malloc((m * m) * sizeof(*colno));
    rowno = (int *) malloc((m * m) * sizeof(*rowno));
    row = (REAL *) malloc((m * m) * sizeof(*row));
    row2 = (REAL *) malloc((2*m * 1) * sizeof(*row2));



    set_add_rowmode(lp, TRUE);

    REAL infinite = get_infinite(lp);
    for (int j = 0; j < m*m; j++) {
        colno[j] = j + 1;
        row[j] = f(j);
        set_bounds(lp, j + 1, 0.0, infinite);
    }
    set_obj_fnex(lp, m*m, row, colno);

    //for (int i = 0; i < 2*m; ++i) {
        for(int j=0; j<m*m; j++) {
            colno[j] = j + 1;
            row[j] = 1.0; //Aeq(i,j);
        }
        add_constraintex(lp, m*m, row, colno, EQ, 1.0); //beq(i)
    //}


    //delete_lp(lp);
    //lprec *lpcopy;// = copy_lp(lp);
    //lpcopy = copy_lp(lp);
    /*for (int i = 0; i < 2*m; ++i) {
        for(int j=0; j<m*m; j++) {
            colno[j] = j + 1;
            row[j] = AA(i,j); //Aeq(i,j);
        }
        //std::cout<<BB(i,0)<<" ";
        add_constraintex(lp, m*m, row, colno, LE, BB(i,0));
    }

    int K = BB.cols();
    Rcpp::NumericVector res(K);
    for (int k = 0; k < K; ++k) {
        lpcopy = copy_lp(lp);
        //set_add_rowmode(lpcopy, F);

        for (int i = 0; i < 2 * m; ++i) {
            for (int j = 0; j < m * m; j++) {
                colno[j] = j + 1;
                row[j] = AA(i, j); //Aeq(i,j);
            }
            //std::cout<<BB(i,0)<<" ";
            add_constraintex(lpcopy, m * m, row, colno, LE, BB(i, k));
        }
        //std::cout<<"\n";
        //set_add_rowmode(lp, FALSE);
        set_add_rowmode(lpcopy, FALSE);


        //set_add_rowmode(lp, FALSE);

        //set_minim(lpcopy);
        //set_verbose(lpcopy, NEUTRAL);

        //if (solve(lp) != OPTIMAL) std::cout << "OPA" << std::endl;

        //
        //res(0) = double(get_objective(lp));
        //std::cout << "obj val = " << double(get_objective(lp)) << std::endl;

        set_minim(lpcopy);
        set_verbose(lpcopy, NEUTRAL);

        if (solve(lpcopy) != OPTIMAL) std::cout << "OPA" << std::endl;
        res(k) = double(get_objective(lpcopy));
        //std::cout << "obj val = " << double(get_objective(lpcopy)) << std::endl;
    }

    return res;

}*/
