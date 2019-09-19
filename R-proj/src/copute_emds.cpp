

#include <Rcpp.h>
#include <RcppEigen.h>
#undef Realloc
#undef Free
#include "lp_lib.h"


// [[Rcpp::export]]
Rcpp::NumericVector emd_Bmat (Rcpp::NumericMatrix B, Rcpp::NumericMatrix A, int m,  Rcpp::NumericVector fin){

    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    MT BB = Rcpp::as<MT>(B);
    //std::cout<<B<<std::endl;

    MT A1(m, m*m);// = matrix(0, m, m * n)
    MT A2(m, m * m);// = matrix(0, n, m * n)

    int k;
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            k = j + i * m;
            A1(i,k) = 1;
            A2(j,k) = 1;
        }
    }

    //MT A(2 * m, m * m);
    //for (int i = 0; i < m; ++i) A.row(i) = A1.row(i);
    //for (int i = m; i < 2*m; ++i) A.row(i) = A2.row(i-m);

    MT AA = Rcpp::as<MT>(A);

    MT Aeq = MT::Ones(2 * m, m * m);
    VT beq = VT::Ones(2 * m);
    VT l = VT::Zero(m * m);
    VT f = Rcpp::as<VT>(fin);
    //std::cout<<AA.rows()<<A.cols()<<std::endl;


    int *colno = NULL, *rowno = NULL;
    REAL *row = NULL, *row2 = NULL;

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

    lpcopy = copy_lp(lp);
    //delete_lp(lp);
    //lprec *lpcopy;// = copy_lp(lp);
    //lpcopy = copy_lp(lp);
    for (int i = 0; i < 2*m; ++i) {
        for(int j=0; j<m*m; j++) {
            colno[j] = j + 1;
            row[j] = AA(i,j); //Aeq(i,j);
        }
        //std::cout<<BB(i,0)<<" ";
        add_constraintex(lp, m*m, row, colno, LE, BB(i,0));
    }

    for (int i = 0; i < 2*m; ++i) {
        for(int j=0; j<m*m; j++) {
            colno[j] = j + 1;
            row[j] = AA(i,j); //Aeq(i,j);
        }
        //std::cout<<BB(i,0)<<" ";
        add_constraintex(lpcopy, m*m, row, colno, LE, BB(i, 1));
    }
    //std::cout<<"\n";
    set_add_rowmode(lp, FALSE);
    set_add_rowmode(lpcopy, FALSE);


    //set_add_rowmode(lp, FALSE);

    set_minim(lp);
    set_verbose(lp, NEUTRAL);

    if (solve(lp) != OPTIMAL) std::cout<<"OPA"<<std::endl;

    Rcpp::NumericVector res(2);
    res(0) = double(get_objective(lp));
    std::cout<<"obj val = "<<double(get_objective(lp))<<std::endl;

    set_minim(lpcopy);
    set_verbose(lpcopy, NEUTRAL);

    if (solve(lpcopy) != OPTIMAL) std::cout<<"OPA"<<std::endl;
    res(1) = double(get_objective(lpcopy));
    std::cout<<"obj val = "<<double(get_objective(lpcopy))<<std::endl;

    return res;

}
