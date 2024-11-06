//
// Created by gregj on 4/30/2024.
//

#ifndef OPTIMATRIXTESTFIXTURE_H
#define OPTIMATRIXTESTFIXTURE_H
#include <memory>

#include "TestFixture.h"
#include "../Adapters/OptimatrixAdapter.h"
#include <vector>

class OptimatrixAdapterTestFixture : public TestFixture{
public:
    bool TestOptimatrixReturnsNotNullValues(const SparseDistanceMatrix* sparse, const ListVector* listVector);
    bool TestOptimatrixClosenessReturnsCorrectValue(const SparseDistanceMatrix* sparse, const ListVector* listVector,
        int expectedSizeOfList);
    bool TestOptimatrixSingletonReturnsCorrectValue(const SparseDistanceMatrix* sparse, const ListVector* listVector,
        int expectedSizeOfList);
    bool TestOptimatrixNameListReturnsCorrectValue(const SparseDistanceMatrix* sparse, const ListVector* listVector,
        int expectedSizeOfList);

private:
    void Setup() override;
    void TearDown() override;
    OptimatrixAdapter* adapter = nullptr;
};



#endif //OPTIMATRIXTESTFIXTURE_H
