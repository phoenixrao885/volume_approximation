// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#include "Eigen/Eigen"
#include "use_double.h"
#include "volume.h"

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

NT factorial(NT n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// Approximating the volume of a convex polytope or body 
// can also be used for integration of concave functions.
// The user should provide the appropriate membership 
// oracles.

int main(const int argc, const char** argv)
{
	//Deafault values
    int n, nexp=1, n_threads=1, W;
    int walk_len,N;
    NT e=1;
    NT exactvol(-1.0);
    bool verbose=false, 
	 rand_only=false, 
	 round_only=false,
	 file=false, 
	 round=false, 
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
         birk=false,
         rotate=false,
         ball_walk=false,
         ball_rad=false,
         experiments=true,
         annealing = false,
         Vpoly=false,
         coordinate=true;
	
	//this is our polytope
	HPolytope<NT> P;
	VPolytope<NT> VP;

	// parameters of CV algorithm
	bool user_W=false, user_N=false, user_ratio=false;
	NT ball_radius=0.0;
	NT C=2.0,ratio,frac=0.1,delta=-1.0,error=0.2;
	
  if(argc<2){
    std::cout<<"Use -h for help"<<std::endl;
    exit(-2);
  }
  
  //parse command line input vars
  for(int i=1;i<argc;++i){
      bool correct=false;

      if(!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")){
          std::cerr<<
                      "Usage:\n"<<
                      "-v, --verbose \n"<<
                      "-rdhr : use random directions HnR, default is coordinate directions HnR\n"
                      "-rand, --rand_only : generates only random points\n"<<
                      "-f1, --file1 [filename_type_Ax<=b] [epsilon] [walk_length] [threads] [num_of_experiments]\n"<<
                      //"-f2, --file2 [filename_type_Ax=b,x>=0] [epsilon] [walk_length] [threads] [num_of_experiments]\n"<<
                      "-fle, --filele : counting linear extensions of a poset\n"<<
                      //"-c, --cube [dimension] [epsilon] [walk length] [threads] [num_of_experiments]\n"<<
                      "--exact : the exact volume\n"<<
                      "--cube : input polytope is a cube\n"<<
                      "-r, --round : enables rounding of the polytope as a preprocess\n"<<
                      "-ro, --round_only : does only rounding to the polytope\n"<<
                      "-e, --error epsilon : the goal error of approximation\n"<<
                      "-w, --walk_len [walk_len] : the random walk length (default 10)\n"<<
                      "-exp [#exps] : number of experiments (default 1)\n"<<
                      "-t, --threads #threads : the number of threads to be used\n"<<
                      "-ΝΝ : use Nearest Neighbor search to compute the boundary oracles\n"<<
                      "-birk_sym : use symmetry to compute more random points (only for Birkhoff polytopes)\n"<<
                      "\n-g_an : use the practical CV algo\n"<<
                      "-w, --walk_len [walk_len] : the random walk length (default 1)\n"<<
                      "-rdhr : use random directions HnR, default is coordinate directions HnR\n"
                      "-e, --error epsilon : the goal error of approximation\n"<<
                      "-bw : use ball walk for sampling\n"<<
                      "-bwr : the radius of the ball walk (default r*chebychev_radius/sqrt(max(1.0, a_i)*dimension\n"<<
                      "-Win : the size of the open window for the ratios convergence\n"<<
                      "-C : a constant for the upper boud of variance/mean^2 in schedule annealing\n"
                      "-N : the number of points to sample in each step of schedule annealing. Default value N = 500*C + dimension^2/2\n"<<
                      "-frac : the fraction of the total error to spend in the first gaussian (default frac=0.1)\n"<<
                      "-ratio : parameter of schedule annealing, larger ratio means larger steps in schedule annealing (default 1-1/dimension)\n"<<
                      std::endl;
          return 0;
      }
      if(!strcmp(argv[i],"--cube")){
          exactvol = std::pow(2,n);
          correct=true;
      }
      if(!strcmp(argv[i],"--exact")){
          exactvol = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-v")||!strcmp(argv[i],"--verbose")){
          verbose=true;
          std::cout<<"Verbose mode\n";
          correct=true;
      }
      if(!strcmp(argv[i],"-rand")||!strcmp(argv[i],"--rand_only")){
          rand_only=true;
          std::cout<<"Generate random points only\n";
          correct=true;
      }
      if(!strcmp(argv[i],"-rdhr")){
          coordinate=false;
          correct=true;
      }
      if(!strcmp(argv[i],"-bw")){
          ball_walk=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-bwr")){
          delta = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-Win")){
          W = atof(argv[++i]);
          user_W=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-ratio")){
          ratio = atof(argv[++i]);
          user_ratio=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-frac")){
          frac = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-C")){
          C = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-N_an")){
          N = atof(argv[++i]);
          user_N=true;
          correct=true;
      }
      //reading from file
      if(!strcmp(argv[i],"-f1")||!strcmp(argv[i],"--file1")){
          file=true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          std::vector<std::vector<NT> > Pin;
          inp.open(argv[++i],std::ifstream::in);
          read_pointset(inp,Pin);
          n = Pin[0][1]-1;
          P.init(Pin);
          if (verbose && P.num_of_hyperplanes()<100){
              std::cout<<"Input polytope: "<<n<<std::endl;
              P.print();
          }
          correct=true;
      }
      if(!strcmp(argv[i],"-f2")||!strcmp(argv[i],"--file2")){
          file=true;
          Vpoly=true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          std::vector<std::vector<NT> > Pin;
          inp.open(argv[++i],std::ifstream::in);
          read_pointset(inp,Pin);
          //std::cout<<"d="<<Pin[0][1]<<std::endl;
          n = Pin[0][1]-1;
          VP.init(Pin);
          if (verbose && VP.num_of_vertices()<100){
              std::cout<<"Input polytope: "<<n<<std::endl;
              VP.print();
          }
          correct=true;
      }
      /*
    if(!strcmp(argv[i],"-f2")||!strcmp(argv[i],"--file2")){
            file=true;
      std::ifstream inp;
      std::vector<std::vector<double> > Pin;
      inp.open(argv[++i],std::ifstream::in);
      read_pointset(inp,Pin);
      //std::cout<<"d="<<Pin[0][1]<<std::endl;
      //n = Pin[0][1]-1;
      P.init(Pin);
      P.rref();
      n=P.dimension();
      //if (verbose && P.num_of_hyperplanes()<1000){
    // std::cout<<"Input polytope: "<<n<<std::endl;
      //  P.print();
      //}
      correct=true;
    }
*/
      //reading linear extensions and order polytopes
      if(!strcmp(argv[i],"-fle")||!strcmp(argv[i],"--filele")){
          file=true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          inp.open(argv[++i],std::ifstream::in);
          std::ofstream os ("order_polytope.ine",std::ofstream::out);
          linear_extensions_to_order_polytope(inp,os);

          std::ifstream inp2;
          inp2.open("order_polytope.ine",std::ifstream::in);
          std::vector<std::vector<NT> > Pin;
          read_pointset(inp2,Pin);
          n = Pin[0][1]-1;
          P.init(Pin);
          std::cout<<"Input polytope: "<<n<<std::endl;
          linear_extensions = true;
          correct=true;
      }
      if(!strcmp(argv[i],"-r")||!strcmp(argv[i],"--round")){
          round = true;
          correct=true;
      }
      if(!strcmp(argv[i],"-e")||!strcmp(argv[i],"--error")){
          e = atof(argv[++i]);
          error=e;
          correct=true;
      }
      if(!strcmp(argv[i],"-w")||!strcmp(argv[i],"--walk_len")){
          walk_len = atof(argv[++i]);
          user_walk_len=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-exp")){
          nexp = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-t")||!strcmp(argv[i],"--threads")){
          n_threads = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-NN")){
          std::cout<<"flann software is needed for this option. Experimental feature."
                  <<"Currently under development."<<std::endl;
          correct=true;
      }
      if(!strcmp(argv[i],"-ro")){
          round_only=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-birk_sym")){
          birk=true;
          correct=true;
      }
      //rotate the polytope randomly
      if(!strcmp(argv[i],"-rot")){
          rotate=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-g_an")){
          annealing=true;
          correct=true;
      }
      if(correct==false){
          std::cerr<<"unknown parameters \'"<<argv[i]<<
                     "\', try "<<argv[0]<<" --help"<<std::endl;
          exit(-2);
      }
      
  }
  
  //Compute chebychev ball//
  std::pair<Point, NT> CheBall;
  double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
  if(!Vpoly) {
      CheBall = P.chebyshev_center();
  }else{
      CheBall = VP.chebyshev_center();
  }
  double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
  if(verbose) std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;
  if(verbose){
      std::cout<<"Chebychev center is: "<<std::endl;
      for(int i=0; i<n; i++){
          std::cout<<CheBall.first[i]<<" ";
      }
      std::cout<<"\nradius is: "<<CheBall.second<<std::endl;
  }
  
  // Set the number of random walk steps
  if(!user_walk_len)
      walk_len = 10 + n / 10;

  if(!user_walk_len) {
      if(!annealing) {
          walk_len = 10 + n / 10;
      }else{
          walk_len = 1;
      }
  }
  if(!user_N)
      N = 500 * ((int) C) + ((int) (n * n / 2));
  if(!user_ratio)
      ratio = 1.0-1.0/(NT(n));
  if(!user_W)
      W = 4*n*n+500;


  // Timings
  double tstart, tstop;

  /* CONSTANTS */
  //error in hit-and-run bisection of P 
  const NT err=0.0000000001;
  const NT err_opt=0.01;

  //bounds for the cube	
  const int lw=0, up=10000, R=up-lw;
  
   /* RANDOM NUMBERS */
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  RNGType rng(seed);
  boost::normal_distribution<> rdist(0,1);
  boost::random::uniform_real_distribution<>(urdist);
  boost::random::uniform_real_distribution<> urdist1(-1,1);

  // If no file specified construct a default polytope
  if(!file){
      P.init(n);
  }

  // If rotate flag is on rotate the polytope
  if(rotate){
      rotating(P);
  }

  // the number of random points to be generated in each K_i
  int rnum = std::pow(e,-2) * 400 * n * std::log(n);
  
  //RUN EXPERIMENTS
  int num_of_exp=nexp;
  double sum_time=0;
  NT min,max,sum=0;
  std::vector<NT> vs;
  NT average, std_dev;
  double Chebtime, sum_Chebtime=double(0);
  NT vol;
  
  for(int i=0; i<num_of_exp; ++i){
      std::cout<<"Experiment "<<i+1<<" ";
      tstart = (double)clock()/(double)CLOCKS_PER_SEC;

      // Setup the parameters
      vars var(rnum,n,walk_len,n_threads,err,e,0,0.0,0,CheBall.second,rng,
               urdist,urdist1,delta,verbose,rand_only,round,NN,birk,ball_walk,coordinate);

      if(round_only){
          // Round the polytope and exit
          std::pair<NT,NT> res_round;
          res_round = rounding_min_ellipsoid(P,CheBall,var);
          NT round_value = res_round.first;
          std::cout<<"\n--------------\nRounded polytope\nH-representation\nbegin\n"<<std::endl;
          P.print();
          std::cout<<"end\n--------------\n"<<std::endl;
      }else{
          // Estimate the volume
          if(annealing){
              // setup the parameters
              vars var2(rnum,n,10 + n/10,n_threads,err,e,0,0.0,0,CheBall.second,rng,
                       urdist,urdist1,delta,verbose,rand_only,round,NN,birk,ball_walk,coordinate);

              vars_g var1(n,walk_len,N,W,1,error,CheBall.second,rng,C,frac,ratio,delta,false,verbose,rand_only,round,NN,birk,ball_walk,coordinate);
              if(!Vpoly) {
                  vol = volume_gaussian_annealing(P, var1, var2, CheBall);
              }else{
                  vol = volume_gaussian_annealing(VP, var1, var2, CheBall);
              }
              tstop = (double)clock()/(double)CLOCKS_PER_SEC;
              std::cout<<"volume computed = "<<vol<<std::endl;
              std::cout<<"Total time = "<<tstop-tstart<<" sec"<<std::endl;
              return 0;
          }
          if(!Vpoly) {
              vol = volume(P, var, var, CheBall);
          }else{
              //vol = volume(VP, var, var, CheBall);
              vol = Vpoly_volume(VP, var);
              NT ext = 1.0/factorial(10);
              std::cout<<"exact volume = "<<ext<<std::endl;
              NT error_vol = std::abs(vol - ext) / ext;
              std::cout<<"error = "<<error_vol<<std::endl;
          }
      }

      NT v1 = vol;

      tstop = (double)clock()/(double)CLOCKS_PER_SEC;

      // Statistics
      sum+=v1;
      if(i==0){max=v1;min=v1;}
      if(v1>max) max=v1;
      if(v1<min) min=v1;
      vs.push_back(v1);
      sum_time +=  tstop-tstart;
      sum_Chebtime += Chebtime;

      if(round)
          std::cout<<" (rounding is ON)";
      std::cout<<std::endl;

      //Compute Statistics
      average=sum/(i+1);
      std_dev=0;
      for(std::vector<NT>::iterator vit=vs.begin(); vit!=vs.end(); ++vit){
          std_dev += std::pow(*vit - average,2);
      }
      std_dev = std::sqrt(std_dev/(i+1));

      std::cout.precision(7);

      //MEMORY USAGE
      //struct proc_t usage;
      //look_up_our_self(&usage);

      //Print statistics
      //std::cout<<"\nSTATISTICS:"<<std::endl;
      if (!experiments){
          std::cout
                  <<"Dimension= "
                  <<n<<" "
                   //<<argv[]<<" "
                  <<"\nNumber of hyperplanes= "
                  <<P.num_of_hyperplanes()<<" "
                  <<"\nNumber of runs= "
                  <<num_of_exp<<" "
                  <<"\nError parameter= "
                  <<e
                  <<"\nTheoretical range of values= "<<" ["
                  <<(1-e)*exactvol<<","
                  <<(1+e)*exactvol<<"] "
                  <<"\nNumber of random points generated in each iteration= "
                  <<rnum<<" "
                  <<"\nRandom walk length= "
                  <<walk_len<<" "
                  <<"\nAverage volume (avg)= "
                  <<average
                  <<"\nmin,max= "
                    " ["
                  <<min<<","
                  <<max<<"] "
                  <<"\nStandard deviation= "
                  <<std_dev<<" "
                  <<"\n(max-min)/avg= "
                  <<(max-min)/average<<" "
                  <<"\nTime(sec)= "
                  <<sum_time/(i+1)<<" "
                  <<"\nTime(sec) Chebyshev= "
                  <<sum_Chebtime/(i+1)<<" "
                    //<<usage.vsize
                  <<std::endl;
    
      if(exactvol!=-1.0){
	      std::cout 
	           <<"\nExact volume= "
	           <<exactvol<<" "
	           <<"\n(vol-avg)/vol= "
	           <<(exactvol-average)/exactvol<<" "
               <<std::endl;
      }
	} else 
    	std::cout 
                 <<n<<" "
                 //<<argv[]<<" "
                 <<P.num_of_hyperplanes()<<" "
                 <<num_of_exp<<" "
                 <<exactvol<<" "
                 <<e<<" ["
                 <<(1-e)*exactvol<<","
                 <<(1+e)*exactvol<<"] "
                 <<rnum<<" "
                 <<walk_len<<" "
                 <<average<<" ["
                 <<min<<","
                 <<max<<"] "
                 <<std_dev<<" "
                 <<(exactvol-average)/exactvol<<" "
                 <<(max-min)/average<<" "
                 <<sum_time/(i+1)<<" "
                 <<sum_Chebtime/(i+1)<<" "
                 //<<usage.vsize
                 <<std::endl;
	}
	
  if(linear_extensions)
		   std::cout <<"Number of linear extensions= "<<vol*factorial(n)<<std::endl;
  
	/*
  // EXACT COMPUTATION WITH POLYMAKE
  /*
	std::ofstream polymakefile;
	polymakefile.open("volume.polymake");
	//print_polymake_volfile(C,polymakefile);
  std::cout<<P[0]<<std::endl;
	print_polymake_volfile2(P,polymakefile);
	system ("polymake volume.polymake");
	std::cout<<std::endl;
  */
  //}
  
  return 0;
}
