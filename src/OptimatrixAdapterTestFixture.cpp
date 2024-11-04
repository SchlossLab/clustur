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
    adapter->ConvertToOptimatrix(xPosition, yPosition, data);
    const bool hasPassed = static_cast<int>(adapter->GetCloseness().size()) == expectedSizeOfList;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixSingletonReturnsCorrectValue(const std::vector<int> &xPosition,
                                                                                const std::vector<int> &yPosition,
                                                                                const std::vector<double> &data,
                                                                                const int expectedSizeOfList) {
    Setup();
    adapter->ConvertToOptimatrix(xPosition, yPosition, data);
    const bool hasPassed = static_cast<int>(adapter->GetSingletons().size()) == expectedSizeOfList;
    TearDown();
    return hasPassed;
}

bool OptimatrixAdapterTestFixture::TestOptimatrixNameListReturnsCorrectValue(const std::vector<int> &xPosition,
                                                                                const std::vector<int> &yPosition,
                                                                                const std::vector<double> &data,
                                                                                const int expectedSizeOfList) {
    Setup();
    adapter->ConvertToOptimatrix(xPosition, yPosition, data);
    const bool hasPassed = static_cast<int>(adapter->GetNameList().size()) == expectedSizeOfList;
    TearDown();
    return hasPassed;
}
bool OptimatrixAdapterTestFixture::ConvertToOptimatrixWithRowData(const std::vector<RowData>& matrixData,
                                                                        const bool sim,
                                                                        const std::vector<std::set<long long>>&
                                                                        expectedResult) {
    Setup();
    auto* optiMatrix = adapter->ConvertToOptimatrix(matrixData, sim);
    //const bool result = expectedResult == optiMatrix->GetCloseness();
    delete optiMatrix;
    TearDown();
    return false;
}

void OptimatrixAdapterTestFixture::Setup() {
    adapter = new OptimatrixAdapter(0.03);
}

void OptimatrixAdapterTestFixture::TearDown() {
    delete(adapter);
}
