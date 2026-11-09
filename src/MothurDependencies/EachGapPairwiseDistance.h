//
// Created by Gregory Johnson on 11/6/26.
//

#ifndef REFACTOR_EACHGAPPAIRWISEDISTANCE_H
#define REFACTOR_EACHGAPPAIRWISEDISTANCE_H
#include "PairwiseDistanceCalculator.h"


class EachGapPairwiseDistance final : public PairwiseDistanceCalculator {
    EachGapPairwiseDistance() = default;
    ~EachGapPairwiseDistance() override = default;
    [[nodiscard]] double Execute(const std::string &sequenceOne, const std::string &sequenceTwo) const override;
};


#endif //REFACTOR_EACHGAPPAIRWISEDISTANCE_H