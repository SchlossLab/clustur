//
// Created by Gregory Johnson on 5/7/24.
//

#include "Tests/OptiDataTestFixture.h"

#include "Adapters/OptimatrixAdapter.h"

bool OptiDataTestFixture::TestGetCloseSequencesReturnCorrectData(const long long index,
                                                                 const std::unordered_set<long long> &expectedResults) {
    Setup();
    const auto resultant = optiData->getCloseSeqs(index);
    TearDown();
    return resultant == expectedResults;
}

bool OptiDataTestFixture::TestIsCloseReturnsIfDataIsClose(const long long index, const long long indexToFind,
    const bool expectedResult) {
    Setup();
    const bool resultant = optiData->isClose(index, indexToFind);
    TearDown();
    return expectedResult == resultant;
}

bool OptiDataTestFixture::TestGetNumCloseReturnsCorrectData(const long long index, const long long expectedResult) {
    Setup();
    const long long resultant = optiData->getNumClose(index);
    TearDown();
    return expectedResult == resultant;
}

bool OptiDataTestFixture::TestGetNameReturnsTheProperName(const long long index, const std::string &expectedResult) {
    Setup();
    const auto resultant = optiData->getName(index);
    TearDown();
    return resultant == expectedResult;
}

bool OptiDataTestFixture::TestGetNumSeqsReturnsCorrectData(const long long expectedResult) {
    Setup();
    const auto resultant = optiData->getNumSeqs();
    TearDown();
    return expectedResult == resultant;
}

bool OptiDataTestFixture::TestGetNumSingletonsReturnsCorrectData(const long long expectedResult) {
    Setup();
    const auto resultant = optiData->getNumSingletons();
    TearDown();
    return expectedResult == resultant;
}

bool OptiDataTestFixture::TestGetListSingleReturnsCorrectData(const int expectedResult) {
    Setup();
    const auto resultant = optiData->getListSingle();
    const auto val = resultant->getNumSeqs();
    TearDown();
    return expectedResult == val;
}

bool OptiDataTestFixture::TestGetCloseRefSeqsReturnsCorrectData(const long long index,
    const std::unordered_set<long long>& expectedResult) {

    Setup();
    const auto resultant = optiData->getCloseRefSeqs(index);

    TearDown();
    return expectedResult == resultant;
}

bool OptiDataTestFixture::TestIsCloseFitReturnsCorrectData(const long long indexOne, const long long indexTwo,
    bool expected, const bool expectedResult) {

    Setup();
    const auto resultant = optiData->isCloseFit(indexOne, indexTwo, expected);
    TearDown();
    return expectedResult == resultant;
}

void OptiDataTestFixture::Setup() {
    const auto xVals = std::vector<int>{0,0,0,1,1,2,3};
    const auto yVals = std::vector<int>{1,2,4,2,4,4,4};
    const auto data = std::vector<float>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
    SparseDistanceMatrix sparse;
    ListVector vec;
    vec.push_back("0");
    vec.push_back("1");
    vec.push_back("2");
    vec.push_back("3");
    vec.push_back("4");
    const size_t size = xVals.size();
    sparse.resize(size);
    for(size_t i = 0; i < size; i++) {
        sparse.addCell(xVals[i], PDistCell(yVals[i], data[i]));
    }
    OptimatrixAdapter adapter(0.03);
    optiData = adapter.ConvertToOptimatrix(&sparse, &vec, false);
}

void OptiDataTestFixture::TearDown() {
    delete(optiData);
}
