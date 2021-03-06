// 
//typedef double                  NT2;

typedef boost::mt19937 RNGType;
typedef CGAL::Gmpq                  NT2;
//typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel2;
//typedef Kernel2::RT					RT;

//typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
//typedef CGAL::Gmpq 			EXACT_NT;
//typedef double NT;
typedef CGAL::Cartesian_d<double> 	      Kernel; 
typedef Kernel::Point_d Point_d;
typedef CGAL::Timer				  timer;
typedef Kernel::Hyperplane_d 		Plane_d;

using namespace Eigen;





NT2 get_numerator(int i, int j, int dim, std::vector<NT2> plane, NT2 z, std::vector<NT2> c){
	
	NT2 q=c[dim], numer, xi, xj, numer2,num;
	int counter, k,m;
	//std::cout<<"[numer]z is: "<<z<<std::endl;
	//std::cout<<"i and j are: "<<i<<" "<<j<<std::endl;
	//MatrixXd A(dim,dim);
	//MatrixXd A2(dim,dim);
	//VectorXd cr(dim);
	//VectorXd x(dim);
	//VectorXd c2(dim);
	
	/*cr(dim-1)=z;
	for (k=0; k<dim; k++){
		c2(k)=c[k];
		A(dim-1,k)=plane[k];
	}
	for (m=0; m<(dim-1); m++){
		for (k=0; k<dim; k++){
			A(m,k)=0.0;
		}
	}*/
	
	if (i==0 || j==0){
		if (i==0){
			
			xj=z/plane[j-1];
			//std::cout<<xj<<std::endl;
			numer=c[j-1]*xj+q;
			num=numer;
			for (k=1; k<dim; k++){
				numer=numer*num;
			}
			//numer=std::pow(numer, ((NT2)dim));
		}else{
			
			
			xi=z/plane[i-1];
			//std::cout<<xi<<std::endl;
			numer=c[i-1]*xi+q;
			num=numer;
			for (k=1; k<dim; k++){
				numer=numer*num;
			}
			//numer=std::pow(numer, ((NT2)dim));
		}
	}else{
		
		
		xi=(z-plane[j-1])/(plane[i-1]-plane[j-1]);
		xj=(z-plane[i-1])/(plane[j-1]-plane[i-1]);
		//std::cout<<xi<<" "<<xj<<std::endl;
		numer=c[i-1]*xi+c[j-1]*xj+q;
		num=numer;
		for (k=1; k<dim; k++){
			numer*=num;
		}
		
	}
	
	return numer;
}






NT2 get_denominator(int i, int j, int dim, std::vector<NT2> plane, std::vector<NT2> c, bool sect){
	
	//std::cout<<"arxh "<<i<<" "<<j<<std::endl;
	int counter,k,m;
	NT2 denom=NT2(1.0),denom2=NT2(1.0),g1,g2;
	std::vector<NT2> g(dim);
	
	
	if (sect){
		
		if (i==0 || j==0){
			//std::cout<<"hello mate!"<<std::endl;
			if (i==0){
				
			
				if (plane[j-1]>0){
					denom=plane[j-1];
				}else{
					denom=-plane[j-1];
				}
				//std::cout<<"det manual is: "<<std::endl;
				//std::cout<<denom<<std::endl;
				g[dim-1]=c[j-1]/plane[j-1];
				//std::cout<<"g[dim-1] is: "<<g[dim-1]<<std::endl;
				denom*=g[dim-1];
				counter=0;
				for (k=0; k<dim; k++){
					if(k!=(j-1)){
						//g[counter]=c[k]-plane[k]*g[dim-1];
						g[counter]=-(c[k]-plane[k]*g[dim-1]);
						//std::cout<<g[counter]<<std::endl;
						denom*=g[counter];
						counter++;
					}
				}
			
			}else{
				
			
				if (plane[i-1]>0){
					denom=plane[i-1];
				}else{
					denom=-plane[i-1];
				}
				//std::cout<<"determinant computed"<<std::endl;
				g[dim-1]=c[i-1]/plane[i-1];
				denom=denom*g[dim-1];
				counter=0;
				for (k=0; k<dim; k++){
					if(k!=(i-1)){
						//g[counter]=c[k]-plane[k]*g[dim-1];
						g[counter]=-(c[k]-plane[k]*g[dim-1]);
						denom*=g[counter];
						counter++;
					}
				}
				
			}
		}else{
			
			
			if(plane[i-1]-plane[j-1]>0){
				denom=(plane[i-1]-plane[j-1]);
			}else{
				denom=(-plane[i-1]+plane[j-1]);
			}
			//denom*=std::abs(plane[i-1]-plane[j-1]);
		//	std::cout<<denom<<std::endl;
			g[dim-2]=(c[i-1]*plane[j-1]-plane[i-1]*c[j-1])/(plane[j-1]-plane[i-1]);
			//std::cout<<"g[dim-2] is: "<<g[dim-2]<<std::endl;
			g[dim-1]=(c[j-1]-c[i-1])/((plane[j-1]-plane[i-1]));
			//std::cout<<"g[dim-1] is: "<<g[dim-1]<<std::endl;
			denom=denom*g[dim-2]*g[dim-1];
			//std::cout<<"denom is: "<<denom<<std::endl;
			counter=0;
			for (k=0; k<dim; k++){
				if (k!=(i-1) && k!=(j-1)){
					//g[counter]=c[k]-g[dim-2]-plane[k]*g[dim-1];
					g[counter]=-(c[k]-g[dim-2]-plane[k]*g[dim-1]);
					//std::cout<<"g[counter] is: "<<g[counter]<<std::endl;
					denom=denom*g[counter];
					counter++;
				}
			}
			
			
		}
		
		
	}else{
		//std::cout<<i<<" "<<j<<std::endl;
		denom=NT2(1.0);
		if (i==0){
			for (k=0; k<dim; k++){
				denom*=-c[k];
			}
		}else{
			g[dim-1]=c[i-1];
			denom=denom*g[dim-1];
			counter=0;
			for (k=0; k<dim; k++){
				if (k!=(i-1)){
					//g[counter]=c[k]-g[dim-1];
					g[counter]=-(c[k]-g[dim-1]);
					denom=denom*g[counter];
					counter++;
				}
			}
		}
		//std::cout<<"denom is: "<<i<<" "<<j<<std::endl;
		//std::cout<<denom<<std::endl;
	}
	
	return denom;
	
}


std::pair< std::vector<int>, std::vector<int> > get_left_right_hyp(int dim, std::vector<NT2> plane, NT2 z){
	
	int i;
	std::vector<int> lefts,rights;
	std::pair< std::vector<int>, std::vector<int> > result;
	
	if (z>0){
		lefts.push_back(0);
	}else{
		rights.push_back(0);
	}
	for (i=0; i<dim; i++){
		if (plane[i]<z){
			lefts.push_back(i+1);
		}else{
			rights.push_back(i+1);
		}
	}
	
	result.first=lefts;
	result.second=rights;
	return result;
	
}







double get_numerator2(int i, int j, int dim, std::vector<double> plane, double z, std::vector<double> c){
	
	double q=c[dim], numer, xi, xj, numer2,num;
	int counter, k,m;
	
	
	if (i==0 || j==0){
		if (i==0){
		
			xj=z/plane[j-1];
			//std::cout<<xj<<std::endl;
			numer=c[j-1]*xj+q;
			num=numer;
			for (k=1; k<dim; k++){
				numer=numer*num;
			}
			//numer=std::pow(numer, ((NT2)dim));
		}else{
			
			xi=z/plane[i-1];
			//std::cout<<xi<<std::endl;
			numer=c[i-1]*xi+q;
			num=numer;
			for (k=1; k<dim; k++){
				numer=numer*num;
			}
			//numer=std::pow(numer, ((NT2)dim));
		}
	}else{
		
		
		xi=(z-plane[j-1])/(plane[i-1]-plane[j-1]);
		xj=(z-plane[i-1])/(plane[j-1]-plane[i-1]);
		//std::cout<<xi<<" "<<xj<<std::endl;
		numer=c[i-1]*xi+c[j-1]*xj+q;
		num=numer;
		for (k=1; k<dim; k++){
			numer*=num;
		}
		
	}
	
	return numer;
}






double get_denominator2(int i, int j, int dim, std::vector<double> plane, std::vector<double> c, bool sect){
	
	//std::cout<<"arxh "<<i<<" "<<j<<std::endl;
	int counter,k,m;
	double denom=1.0,denom2=1.0,g1,g2;
	std::vector<double> g(dim);
	
	
	if (sect){
		
		if (i==0 || j==0){
			//std::cout<<"hello mate!"<<std::endl;
			if (i==0){
				
			
				if (plane[j-1]>0){
					denom=plane[j-1];
				}else{
					denom=-plane[j-1];
				}
				//std::cout<<"det manual is: "<<std::endl;
				//std::cout<<denom<<std::endl;
				g[dim-1]=c[j-1]/plane[j-1];
				//std::cout<<"g[dim-1] is: "<<g[dim-1]<<std::endl;
				denom*=g[dim-1];
				counter=0;
				for (k=0; k<dim; k++){
					if(k!=(j-1)){
						//g[counter]=c[k]-plane[k]*g[dim-1];
						g[counter]=-(c[k]-plane[k]*g[dim-1]);
						//std::cout<<g[counter]<<std::endl;
						denom*=g[counter];
						counter++;
					}
				}
			//	std::cout<<"denom is: "<<i<<" "<<j<<std::endl;
			//	std::cout<<denom<<std::endl;
				//std::cout<<"denom2 is: "<<std::endl;
				//std::cout<<denom2<<std::endl;
			}else{
				
			
				if (plane[i-1]>0){
					denom=plane[i-1];
				}else{
					denom=-plane[i-1];
				}
				//std::cout<<"determinant computed"<<std::endl;
				g[dim-1]=c[i-1]/plane[i-1];
				denom=denom*g[dim-1];
				counter=0;
				for (k=0; k<dim; k++){
					if(k!=(i-1)){
						//g[counter]=c[k]-plane[k]*g[dim-1];
						g[counter]=-(c[k]-plane[k]*g[dim-1]);
						denom*=g[counter];
						counter++;
					}
				}
				//std::cout<<"denom is: "<<i<<" "<<j<<std::endl;
			//	std::cout<<denom<<std::endl;
				//std::cout<<"denom2 is: "<<std::endl;
				//std::cout<<denom2<<std::endl;
			}
		}else{
			
		
			if(plane[i-1]-plane[j-1]>0){
				denom=(plane[i-1]-plane[j-1]);
			}else{
				denom=(-plane[i-1]+plane[j-1]);
			}
			//denom*=std::abs(plane[i-1]-plane[j-1]);
		//	std::cout<<denom<<std::endl;
			g[dim-2]=(c[i-1]*plane[j-1]-plane[i-1]*c[j-1])/(plane[j-1]-plane[i-1]);
			//std::cout<<"g[dim-2] is: "<<g[dim-2]<<std::endl;
			g[dim-1]=(c[j-1]-c[i-1])/((plane[j-1]-plane[i-1]));
			//std::cout<<"g[dim-1] is: "<<g[dim-1]<<std::endl;
			denom=denom*g[dim-2]*g[dim-1];
			//std::cout<<"denom is: "<<denom<<std::endl;
			counter=0;
			for (k=0; k<dim; k++){
				if (k!=(i-1) && k!=(j-1)){
					//g[counter]=c[k]-g[dim-2]-plane[k]*g[dim-1];
					g[counter]=-(c[k]-g[dim-2]-plane[k]*g[dim-1]);
					//std::cout<<"g[counter] is: "<<g[counter]<<std::endl;
					denom=denom*g[counter];
					counter++;
				}
			}
			//std::cout<<"denom is: "<<i<<" "<<j<<std::endl;
			//std::cout<<denom<<std::endl;
			//std::cout<<"denom2 is: "<<std::endl;
			//std::cout<<denom2<<std::endl;
			
		}
		
		
	}else{
		//std::cout<<i<<" "<<j<<std::endl;
		denom=1.0;
		if (i==0){
			for (k=0; k<dim; k++){
				denom*=-c[k];
			}
		}else{
			g[dim-1]=c[i-1];
			denom=denom*g[dim-1];
			counter=0;
			for (k=0; k<dim; k++){
				if (k!=(i-1)){
					//g[counter]=c[k]-g[dim-1];
					g[counter]=-(c[k]-g[dim-1]);
					denom=denom*g[counter];
					counter++;
				}
			}
		}
		//std::cout<<"denom is: "<<i<<" "<<j<<std::endl;
		//std::cout<<denom<<std::endl;
	}
	
	return denom;
	
}


std::pair< std::vector<int>, std::vector<int> > get_left_right_hyp2(int dim, std::vector<double> plane, double z){
	
	int i;
	std::vector<int> lefts,rights;
	std::pair< std::vector<int>, std::vector<int> > result;
	
	if (z>0){
		lefts.push_back(0);
	}else{
		rights.push_back(0);
	}
	for (i=0; i<dim; i++){
		if (plane[i]<z){
			lefts.push_back(i+1);
		}else{
			rights.push_back(i+1);
		}
	}
	
	result.first=lefts;
	result.second=rights;
	return result;
	
}
