//
// Created by Gregory Johnson on 8/13/26.
//

#ifndef REFACTOR_CREATEDISTANCEMATRIX_H
#define REFACTOR_CREATEDISTANCEMATRIX_H
#include "PairwiseDistanceCalculator.h"
#include "../DataStructures/SparseDistanceMatrix.h"


class CreateDistanceMatrix {
public:
    static SparseDistanceMatrix CreateSparseDistanceMatrix(const std::vector<std::string>& sequences,
        PairwiseDistanceCalculator* calculator, double cutoff);
};


#endif //REFACTOR_CREATEDISTANCEMATRIX_H