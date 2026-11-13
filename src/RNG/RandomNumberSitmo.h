//
// Created by Gregory Johnson on 9/18/26.
//

#ifndef RANDOMNUMBERSITMO_H
#define RANDOMNUMBERSITMO_H
#include <threefry.h>
class RandomNumberSitmo {
public:
    using result_type = long long;
    explicit RandomNumberSitmo(int seed) {
        rngEngine.seed(seed);
    }
    static constexpr result_type min()
    {
        return 0;
    }

    static constexpr result_type max()
    {
        return RAND_MAX;
    }

    result_type operator()()
    {
        return NextRandomValue();
    }
    // [[Rcpp::depends(sitmo)]]
    result_type NextRandomValue() {
        constexpr double dis = max() - min();
        const double result = min() + (static_cast<double>(rngEngine()) / (sitmo::threefry::max())) * (dis);
        return static_cast<long long>(result);
    }
    result_type GetRandomNumber(const double min = 0.0, const double max = 1.0) {
        const double dis = max - min + 1;
        return min + (static_cast<double>(rngEngine()) / (sitmo::threefry::max())) * (dis);
    }

private:
    sitmo::threefry rngEngine;
};
#endif //RANDOMNUMBERSITMO_H