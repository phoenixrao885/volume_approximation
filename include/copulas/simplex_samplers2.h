// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

#ifndef SIMPLEX_SAMPLERS2_H
#define SIMPLEX_SAMPLERS2_H

//#include <boost/random/uniform_real_distribution.hpp>

template <typename NT, class RNGType, class MT>
void exp_simplex (int d, int N, MT &points) {

    boost::random::uniform_real_distribution<> urdist(0, 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    NT sum;

    for (int i = 0; i < N; ++i) {
        sum = 0.0;
        for (int j = 0; j < d; ++j) {
            points(j,i) = -log(urdist(rng));
            sum += points(j,i);
        }
        points.col(i) *= 1.0 / sum;
    }
}

#endif
