//
// Created by Gregory Johnson on 9/18/26.
//

#ifndef RANDOMNUMBERSITMO_H
#define RANDOMNUMBERSITMO_H
#include <threefry.h>
class RandomNumberSitmo {
public:
    using result_type = uint64_t;
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
        double dis = max() - min();
        double result = min() + (static_cast<double>(rngEngine()) / (sitmo::threefry::max())) * (dis);
        return static_cast<uint64_t>(result);
    }
    result_type GetRandomNumber(const double min = 0.0, const double max = 1.0) {
        double dis = max - min + 1;
        return min + (static_cast<double>(rngEngine()) / (sitmo::threefry::max())) * (dis);
    }

private:
    sitmo::threefry rngEngine;
};
#endif //RANDOMNUMBERSITMO_H