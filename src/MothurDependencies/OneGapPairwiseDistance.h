//
// Created by Gregory Johnson on 8/11/26.
//

#ifndef REFACTOR_ONEGAPPAIRWISEDISTANCE_H
#define REFACTOR_ONEGAPPAIRWISEDISTANCE_H
#include "PairwiseDistanceCalculator.h"


class OneGapPairwiseDistance final : public PairwiseDistanceCalculator {
public:
    OneGapPairwiseDistance() = default;
    ~OneGapPairwiseDistance() override = default;
    double Execute(const std::string &sequenceOne, const std::string &sequenceTwo) override;
};


#endif //REFACTOR_ONEGAPPAIRWISEDISTANCE_H