//
// Created by gregj on 4/30/2024.
//

#include "Tests/OptimatrixAdapterTestFixture.h"


bool OptimatrixAdapterTestFixture::TestOptimatrixReturnsCorrectValue(const std::vector<int> &xPosition,
                                                                     const std::vector<int> &yPosition,
                                                                     const std::vector<double> &data) {
    Setup();
    const bool hasPassed = adapter->ConvertToOptimatrix(xPosition, yPosition, data) != nullptr;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixClosenessReturnsCorrectValue(const std::vector<int> &xPosition,
                                                                              const std::vector<int> &yPosition,
                                                                              const std::vector<double> &data,
                                                                              const int expectedSizeOfList) {
    Setup();
    bool hasPassed = false;
    adapter->ConvertToOptimatrix(xPosition, yPosition, data);
    hasPassed = adapter->GetCloseness().size() == expectedSizeOfList;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixSingletonReturnsCorrectValue(const std::vector<int> &xPosition,
                                                                                const std::vector<int> &yPosition,
                                                                                const std::vector<double> &data,
                                                                                const int expectedSizeOfList) {
    Setup();
    bool hasPassed = false;
    adapter->ConvertToOptimatrix(xPosition, yPosition, data);
    hasPassed = adapter->GetSingletons().size() == expectedSizeOfList;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixNameListReturnsCorrectValue(const std::vector<int> &xPosition,
                                                                                const std::vector<int> &yPosition,
                                                                                const std::vector<double> &data,
                                                                                const int expectedSizeOfList) {

    // OptimatrixAdapter adapt(0.03);
    // adapt.ConvertToOptimatrix(xVals, yVals, data);
    Setup();
    bool hasPassed = false;
    adapter->ConvertToOptimatrix(xPosition, yPosition, data);
    hasPassed = adapter->GetNameList().size() == expectedSizeOfList;
    TearDown();
    return hasPassed;
}

void OptimatrixAdapterTestFixture::Setup() {
    adapter = std::make_unique<OptimatrixAdapter>(0.03);
}

void OptimatrixAdapterTestFixture::TearDown() {
    adapter.reset();
}
