//
// Created by Gregory Johnson on 11/6/26.
//

#ifndef REFACTOR_EACHGAPPAIRWISEDISTANCE_H
#define REFACTOR_EACHGAPPAIRWISEDISTANCE_H
#include "PairwiseDistanceCalculator.h"


class EachGapPairwiseDistance final : public PairwiseDistanceCalculator {
public:
    EachGapPairwiseDistance() = default;
    ~EachGapPairwiseDistance() override = default;
    [[nodiscard]] double Execute(const std::string &sequenceOne, const std::string &sequenceTwo,
        double cutoff = 1) const override;
};


#endif //REFACTOR_EACHGAPPAIRWISEDISTANCE_H