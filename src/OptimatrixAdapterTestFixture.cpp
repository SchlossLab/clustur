//
// Created by gregj on 4/30/2024.
//

#include "Tests/OptimatrixAdapterTestFixture.h"


bool OptimatrixAdapterTestFixture::TestOptimatrixReturnsNotNullValues(const SparseDistanceMatrix* sparse,
                                                                    const ListVector* listVector) {
    Setup();
    const bool hasPassed = adapter->ConvertToOptimatrix(sparse, listVector, false) != nullptr;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixClosenessReturnsCorrectValue(const SparseDistanceMatrix* sparse,
                                                                              const ListVector* listVector,
                                                                              const int expectedSizeOfList) {
    Setup();
    const auto optiMatrix = adapter->ConvertToOptimatrix(sparse, listVector, false);
    const bool hasPassed = static_cast<int>(optiMatrix->GetCloseness().size()) == expectedSizeOfList;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixSingletonReturnsCorrectValue(const SparseDistanceMatrix* sparse,
                                                                              const ListVector* listVector,
                                                                              const int expectedSizeOfList) {
    Setup();
    const auto optiMatrix = adapter->ConvertToOptimatrix(sparse, listVector, false);
    const bool hasPassed = static_cast<int>(optiMatrix->GetSingletons().size()) == expectedSizeOfList;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixNameListReturnsCorrectValue(const SparseDistanceMatrix* sparse,
                                                                              const ListVector* listVector,
                                                                              const int expectedSizeOfList) {
    Setup();
    const auto optiMatrix = adapter->ConvertToOptimatrix(sparse, listVector, false);
    const bool hasPassed = static_cast<int>(optiMatrix->GetNameList().size()) == expectedSizeOfList;
    TearDown();
    return hasPassed;
}
void OptimatrixAdapterTestFixture::Setup() {
    adapter = new OptimatrixAdapter(0.03);
}

void OptimatrixAdapterTestFixture::TearDown() {
    delete(adapter);
}
