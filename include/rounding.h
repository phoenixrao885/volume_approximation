// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ROUNDING_H
#define ROUNDING_H

//ROUNDING
template <class MT>
MT getPointsMat(std::list<Point> randPoints, int dim){

    MT S(randPoints.size(),dim);
    for(unsigned int i=0; i<randPoints.size(); i++){
        Point p=randPoints.front();
        randPoints.pop_front();
        for (unsigned int j=0; j<dim; j++){
            S(i,j)=p[j];
        }
    }
    
    return S;
}


template <class T1>
std::pair<Point, NT> approx_R(T1 &P, vars var){
    std::pair<Point,double> InnerBall=P.ComputeInnerBall();
    Point c=InnerBall.first;
    NT radius = InnerBall.second;

    int n=var.n, walk_len=var.walk_steps;
    //Random_points_on_sphere_d<Point> gen (n, radius);
    //Point p = gen.sample_point(var.rng);
    Point p = get_point_on_Dsphere(n, radius);
    p = p + c;
    std::list<Point> randPoints; //ds for storing rand points
    //use a large walk length e.g. 1000
    rand_point_generator(P, p, 1, 50*n, randPoints, var);
    //if (print) std::cout<<"First random point: "<<p<<std::endl;

    // 3. Sample points from P
    //randPoints.push_front(p);
    int num_of_samples = std::pow(1.0,-2) * 400 * n * std::log(n);;//this is the number of sample points will used to compute min_ellipoid
    //if(print) std::cout<<"\nCompute "<<num_of_samples<<" random points in P"<<std::endl;
    rand_point_generator(P, p, num_of_samples*10, 1, randPoints, var);
    NT current_dist, max_dist;
    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist){
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    NT R=max_dist/radius;
    return std::pair<Point,NT> (c,R);
}


//Needs developing (experimental)
template <class T1>
NT rounding_SVD(T1 &P , Point c, NT radius, vars &var){
    typedef typename T1::MT 	MT;
    typedef typename T1::VT 	VT;
    int n=var.n, walk_len=var.walk_steps;
    bool print=var.verbose;
    // 1. Compute the Chebychev ball (largest inscribed ball) with center and radius 
	//Point c(n);       //center
    //K radius;
    //P.chebyshev_center(c,radius);
    
    // 2. Generate the first random point in P
  // Perform random walk on random point in the Chebychev ball 
	//if (print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
	vars var2=var;
	var2.coordinate=false;
    Point p = get_point_on_Dsphere(n, radius);
	p = p + c;
	std::list<Point> randPoints; //ds for storing rand points
	//use a large walk length e.g. 1000
	rand_point_generator(P, p, 1, 50*n, randPoints, var2);
	//if (print) std::cout<<"First random point: "<<p<<std::endl;
    // 3. Sample points from P
	//randPoints.push_front(p);
	int num_of_samples = std::pow(1.0,-2) * 400 * n * std::log(n);;//this is the number of sample points will used to compute min_ellipoid
	//if(print) std::cout<<"\nCompute "<<num_of_samples<<" random points in P"<<std::endl;
	rand_point_generator(P, p, num_of_samples*10, 1, randPoints, var2);
    NT current_dist, max_dist;
    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist){
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    NT R=max_dist/radius;
    if(print) std::cout<<"R = "<<max_dist<<" r = "<<radius<<"ratio R/r = "<<R<<"\n"<<std::endl;
    
    // 4. Compute the transformation matrix T
    MT T = Eigen::MatrixXd::Identity(n,n);
    bool well_rounded=false;
    int t=8*n*n*n;
    //int t=var.m;
    int tries=0;
    MT S=Eigen::MatrixXd::Identity(n,n);
    std::pair<Point,NT> res;
    while(!well_rounded){
        tries++;
        randPoints.clear();
        T1 P2(P);
        P2.linear_transformIt(T.inverse());   //We have to sample from the transformed body
        res=solveLP(P2.get_matrix(), P2.dimension());
        c=res.first;
        p = get_point_on_Dsphere(n, res.second);
        p = p + c;
        rand_point_generator(P2, p, 1, 50*n, randPoints, var2);
        randPoints.clear();
        rand_point_generator(P2, p, t, 1, randPoints, var2);
        MT PM=getPointsMat<MT>(randPoints,n);
        Eigen::JacobiSVD<MT> svd(PM, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if(print) std::cout<<svd.singularValues()<<"\n"<<std::endl;
        NT min=svd.singularValues()(0);
        for(unsigned int i=1; i<n; i++){
            if(svd.singularValues()(i)<min){
                min=svd.singularValues()(i);
            }
        }
        for(unsigned int i=0; i<n; i++){
            S(i,i)=svd.singularValues()(i)/min;
        }
        if(print) std::cout<<S<<"\n"<<std::endl;

        T=svd.matrixV()*S.inverse()*T;
        if(print) std::cout<<T<<"\n"<<std::endl;
        well_rounded=true;
        for (unsigned int i=0; i<n; i++){
            if (S(i,i)>2.0){
                if (tries>((int)log2(R))){
                    t=t*2;
                    tries=0;
                }
                well_rounded=false;
                break;
            }
        }
        if (well_rounded){
            P.linear_transformIt(T.inverse());
        }
    }
    
    return std::abs(T.determinant());
}


// ----- ROUNDING ------ //
// main rounding function
template <class T1, typename FT>

std::pair <FT, FT> rounding_min_ellipsoid(T1 &P , std::pair<Point,FT> InnerBall, vars &var) {

    typedef typename T1::MT 	MT;
    typedef typename T1::VT 	VT;
    int n=var.n, walk_len=var.walk_steps, i, j = 0;
    bool print=var.verbose;
    Point c = InnerBall.first;
    FT radius = InnerBall.second, R = 1;
    std::list<Point> randPoints; //ds for storing rand points
    if (!P.get_points_for_rounding(randPoints)) {  // If P is a V-polytope then it will store its vertices in randPoints
        // If P is not a V-Polytope or number_of_vertices>20*domension
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        Point p = get_point_on_Dsphere(n, radius);
        p = p + c;

        //use a large walk length e.g. 1000
        rand_point_generator(P, p, 1, 50*n, randPoints, var);
        // 3. Sample points from P
        int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
        randPoints.clear();
        rand_point_generator(P, p, num_of_samples, walk_len, randPoints, var);
        /*FT current_dist, max_dist;
        for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
            current_dist=(*pit-c).squared_length();
            if(current_dist>max_dist){
                max_dist=current_dist;
            }
        }
        max_dist=std::sqrt(max_dist);
        R=max_dist/radius;*/
    }

    // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
    boost::numeric::ublas::matrix<double> Ap(n,randPoints.size());
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<FT>::iterator qit;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            Ap(i,j)=double(*qit);
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n);
    boost::numeric::ublas::vector<double> c2(n);
    size_t w=1000;
    FT elleps=Minim::KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

    MT E(n,n);
    VT e(n);

    //Get ellipsoid matrix and center as Eigen objects
    for(unsigned int i=0; i<n; i++){
        e(i)=FT(c2(i));
        for (unsigned int j=0; j<n; j++){
            E(i,j)=FT(Q(i,j));
        }
    }


    //Find the smallest and the largest axes of the elliposoid
    Eigen::EigenSolver<MT> eigensolver(E);
    FT rel = std::real(eigensolver.eigenvalues()[0]);
    FT Rel = std::real(eigensolver.eigenvalues()[0]);
    for(unsigned int i=1; i<n; i++){
        if(std::real(eigensolver.eigenvalues()[i])<rel) rel=std::real(eigensolver.eigenvalues()[i]);
        if(std::real(eigensolver.eigenvalues()[i])>Rel) Rel=std::real(eigensolver.eigenvalues()[i]);
    }

    Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
    MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    //Shift polytope in order to contain the origin (center of the ellipsoid)
    P.shift(e);

    MT L_1 = L.inverse();
    P.linear_transformIt(L_1.transpose());

    return std::pair<FT,FT> (L_1.determinant(),rel/Rel);
}


template <class MT, class VT, class T1>
std::pair< std::pair<MT, VT>, NT> comp_lintrans(T1 &P, vars &var) {

    int n=var.n, walk_len=var.walk_steps, i, j = 0;
    bool print=var.verbose, get_points = true;
    std::pair<Point, NT> CheBall = P.ComputeInnerBall();
    Point c = CheBall.first;
    NT radius = CheBall.second, ratio;
    std::list<Point> randPoints; //ds for storing rand points
    int num_of_samples;// = 10*n;//this is the number of sample points will used to compute min_ellipoid
    if (!P.get_points_for_rounding(randPoints)) {  // If P is a V-polytope then it will store its vertices in randPoints
        get_points = false; // P is H polytope or V-polytope with num_of_vertices>20*dimension
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        Point p = get_point_on_Dsphere(n, radius);
        p = p + c;
        //use a large walk length e.g. 1000
        rand_point_generator(P, p, 1, 50*n, randPoints, var);
        // 3. Sample points from P
        randPoints.clear();
        rand_point_generator(P, p, 100*n, walk_len, randPoints, var); // sample >10*n in order to estimate ratio between radius of largest inscribed ball and smallest enclosing

        // estimate ratio between radius of largest inscribed ball and smallest enclosing for the stopping criterion
        NT current_dist, max_dist;
        for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
            current_dist=(*pit-c).squared_length();
            if(current_dist>max_dist){
                max_dist=current_dist;
            }
        }
        max_dist=std::sqrt(max_dist);
        ratio=radius/max_dist;

    }
    if (get_points) {
        num_of_samples = randPoints.size(); // P is V-polytope
    } else {
        num_of_samples = 10*n; // use only 10*n for the minimum volume enclosing ellipsoid
    }

    // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
    boost::numeric::ublas::matrix<double> Ap(n,num_of_samples);
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    for ( j=0; j<num_of_samples; rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            Ap(i,j)=double(*qit);
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n);
    boost::numeric::ublas::vector<double> c2(n);
    size_t w=1000;
    NT elleps=Minim::KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

    MT E(n,n);
    VT e(n);

    //Get ellipsoid matrix and center as Eigen objects
    for(int i=0; i<n; i++){
        e(i)=NT(c2(i));
        for (int j=0; j<n; j++){
            E(i,j)=NT(Q(i,j));
        }
    }

    if (get_points) { // V-polytope. Use the ratio between minimum and maximum axe of the enclosing ellipsoid as stopping criterion
        //Find the smallest and the largest axes of the elliposoid
        Eigen::EigenSolver <MT> eigensolver(E);
        NT rel = std::real(eigensolver.eigenvalues()[0]);
        NT Rel = std::real(eigensolver.eigenvalues()[0]);
        for (int i = 1; i < n; i++) {
            if (std::real(eigensolver.eigenvalues()[i]) < rel) rel = std::real(eigensolver.eigenvalues()[i]);
            if (std::real(eigensolver.eigenvalues()[i]) > Rel) Rel = std::real(eigensolver.eigenvalues()[i]);
        }
        ratio = rel/Rel;
    }

    std::pair<MT, VT> Ell;
    Ell.first = E; Ell.second = e;

    return std::pair< std::pair<MT, VT>, NT> (Ell, ratio);
}


template <class T1, class MT, class VT, typename FT>
std::pair< std::pair<MT, VT>, FT> is_last_rounding(T1 &P2,
                                                   std::pair<MT, VT> Ell,
                                                   FT ratio,
                                                   FT &round_value,
                                                   bool &comp_next,
                                                   vars &var) {

    T1 P3(P2);
    MT E = Ell.first;
    VT e = Ell.second;

    Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
    MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    //Shift polytope in order to contain the origin (center of the ellipsoid)
    P3.shift(e);

    MT L_1 = L.inverse();
    // apply linear transformation in test polytope P3 in order to check the new ratio
    P3.linear_transformIt(L_1.transpose());

    std::pair< std::pair<MT, VT>, FT> res = comp_lintrans<MT,VT>(P3, var);

    //std::cout<<"now = "<<ratio<<" next = "<<res.second<<std::endl;

    if (res.second > ratio) {  // if the ratio is largest then accept linear transformation
        P2 = P3;
        round_value *= L_1.determinant();
    } else { // if the ratio is smaller then stop rounding
        comp_next = false;
        return std::pair< std::pair<MT, VT>, FT> (Ell,ratio);
    }

    return res;

}


template <class T1>
std::pair<NT,NT> round_polytope(T1 &P, vars &var) {

    typedef typename T1::MT 	MT;
    typedef typename T1::VT 	VT;
    std::pair< std::pair<MT, VT>, NT> res;
    res = comp_lintrans<MT, VT>(P, var);
    bool comp_next = true;
    NT round_value = NT(1.0);
    T1 P2(P);

    while (comp_next) {
        //std::cout<<"round_value = "<<round_value<<std::endl;
        // comp_next is false when ratio gets smaller for the first time
        res = is_last_rounding(P2, res.first, res.second, round_value, comp_next, var);
        //std::cout<<"round_value = "<<round_value<<"\n"<<std::endl;
    }
    P = P2;
    return std::pair<NT, NT> (round_value, res.second);


}


// -------- ROTATION ---------- //
/*
template <class T>
double rotating_old(T &P){

  bool print = true; 
  if(print) std::cout<<"\nRotate..."<<std::endl;
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MT;
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> VT;
  
  int m=P.num_of_hyperplanes();
  int d=P.dimension();
  
  MT A(m,d);
  VT b(m);
  for(int i=0; i<m; ++i){
		b(i) = P.get_coeff(i,0);
		for(int j=1; j<d+1; ++j){
		  A(i,j-1) = P.get_coeff(i,j);
		}
	}
  std::cout<<A<<"\n"<<b<<std::endl;
  
  // Construct rotation matrix by 30o
  MT R(m,m);
  for(int i=0; i<m; ++i){
		for(int j=0; j<m; ++j){
		  R(i,j) = 0;
		  if (i==j) R(i,j) = 1;
		  if (i==0 && j==0) R(i,j) = std::sqrt(3)/2;
		  if (i==1 && j==0) R(i,j) = -1.0/2.0;
		  if (i==0 && j==1) R(i,j) = 1.0/2.0;
		  if (i==1 && j==1) R(i,j) = std::sqrt(3)/2;
		  if (i==d && j==d) R(i,j) = std::sqrt(3)/2;
		  if (i==d+1 && j==d) R(i,j) = -1.0/2.0;
		  if (i==d && j==d+1) R(i,j) = 1.0/2.0;
		  if (i==d+1 && j==d+1) R(i,j) = std::sqrt(3)/2;
		}
	}
	//std::cout<<R<<std::endl;
	
	A = R*A;
	//b = R*b;
	
	std::cout<<A<<"\n"<<b<<std::endl;
	
	// Write changes (actually perform rotation) to the polytope!
	for(int i=0; i<m; ++i){
		P.put_coeff(i,0,b(i));
		for(int j=1; j<d+1; ++j){
		  P.put_coeff(i,j,A(i,j-1));
		}
	}
	
  std::cout<<R.determinant()<<"\n"<<b<<std::endl;
  
	return R.determinant();
}*/

// -------- ROTATION ---------- //
template <class T>
NT rotating(T &P){

  typedef typename T::MT 	MT;
  typedef typename T::VT 	VT;

  int n = P.dimension();

  // pick a random rotation
  MT R = MT::Random(n,n);
  Eigen::JacobiSVD<MT> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

  // apply rotation to the polytope P
  P.linear_transformIt(svd.matrixU());

  return std::abs(svd.matrixU().inverse().determinant());
}

#endif
