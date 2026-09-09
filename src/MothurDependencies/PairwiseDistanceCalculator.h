//
// Created by Gregory Johnson on 8/11/26.
//

#ifndef REFACTOR_PAIRWISEDISTANCECALCULATOR_H
#define REFACTOR_PAIRWISEDISTANCECALCULATOR_H
#include <string>
#include <vector>

class PairwiseDistanceCalculator {
public:
    PairwiseDistanceCalculator() = default;
    virtual ~PairwiseDistanceCalculator() = default;
    virtual double Execute(const std::string& sequenceOne, const std::string& sequenceTwo) const = 0;
};


#endif //REFACTOR_PAIRWISEDISTANCECALCULATOR_H