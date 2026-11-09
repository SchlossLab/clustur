//
// Created by Gregory Johnson on 11/6/26.
//

#ifndef REFACTOR_ONEGAPIGNORE_H
#define REFACTOR_ONEGAPIGNORE_H
#include "PairwiseDistanceCalculator.h"


class OneGapIgnore final : public PairwiseDistanceCalculator{
public:
    OneGapIgnore() = default;
    ~OneGapIgnore() override = default;
    [[nodiscard]] double Execute(const std::string &sequenceOne, const std::string &sequenceTwo) const override;
};


#endif //REFACTOR_ONEGAPIGNORE_H