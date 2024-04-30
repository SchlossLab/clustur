//
// Created by gregj on 4/30/2024.
//

#include "Tests/OptimatrixAdapterTestFixture.h"


bool OptimatrixAdapterTestFixture::TestOptimatrixReturnsCorrectValue(const std::vector<int> &xPosition,
                                                                     const std::vector<int> &yPosition, const std::vector<double> &data) {
    Setup();
    const bool hasPassed = adapter->ConvertToOptimatrix(xPosition, yPosition, data) != nullptr;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixClosenessReturnsCorrectValue() {
    return false;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixSingletonReturnsCorrectValue() {
    return false;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixNameListReturnsCorrectValue() {
    return false;
}

void OptimatrixAdapterTestFixture::Setup() {
    adapter = std::make_unique<OptimatrixAdapter>(0.2);
}

void OptimatrixAdapterTestFixture::TearDown() {
    adapter.reset();
}
