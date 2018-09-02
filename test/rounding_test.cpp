// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template <typename NT, class RNGType, class Polytope>
void rounding_test(Polytope &P, bool rot, NT expected, NT tolerance=0.1)
{

    typedef typename Polytope::PolytopePoint Point;

    // Setup the parameters
    int n = P.dimension();
    int walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,true);

    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    //apply rotation if requested
    NT rot_val;
    if(rot){
        std::cout << "\n--- Rotation is ON "<< std::endl;
        rot_val = rotating<NT>(P);
        std::cout << "Rotation value = "<<rot_val<<std::endl;
    }

    std::pair<NT,NT> res_round = round_polytope<NT>(P, var);
    NT round_value = res_round.first;
    std::pair<Point,NT> CheBall;
    CheBall = P.ComputeInnerBall();
    std::cout<<"\nround value is: "<<round_value<<std::endl;

    //estimate volume
    NT vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        vol += round_value*volume(P,var,var,CheBall);
    }
    NT error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    CHECK(error < tolerance);

}


template <typename NT>
void call_test_rot_skinny_cubes() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing rounding of rotated H-skinny_cube10" << std::endl;
    P = gen_skinny_cube<Hpolytope>(10);
    rounding_test<NT, RNGType>(P, true, 102400.0);

    //std::cout << "\n--- Testing rounding of rotated H-skinny_cube20" << std::endl;
    //P = gen_skinny_cube<Hpolytope>(20);
    //rounding_test<NT, RNGType>(P, true, 104857600.0, 0.2);

}

template <typename NT>
void call_test_skinny_cubes() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing rounding of H-skinny_cube10" << std::endl;
    P = gen_skinny_cube<Hpolytope>(10);
    rounding_test<NT, RNGType>(P, false, 102400.0);

    std::cout << "\n--- Testing rounding of H-skinny_cube20" << std::endl;
    P = gen_skinny_cube<Hpolytope>(20);
    rounding_test<NT, RNGType>(P, false, 104857600.0);

}


TEST_CASE("round_rot_skinny_cube") {
    call_test_rot_skinny_cubes<double>();
    call_test_rot_skinny_cubes<float>();
    call_test_rot_skinny_cubes<long double>();
}

TEST_CASE("round_skinny_cube") {
    call_test_skinny_cubes<double>();
    call_test_skinny_cubes<float>();
    call_test_skinny_cubes<long double>();
}
