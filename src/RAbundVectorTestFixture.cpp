//
// Created by Gregory Johnson on 7/12/24.
//

#include "Tests/RAbundVectorTestFixture.h"

#include "Adapters/MatrixAdapter.h"

bool RAbundVectorTestFixture::TestRemove(const int bin, const int expectedResult) {
    Setup();
    const int abund = radbund.remove(bin);
    TearDown();
    return abund == expectedResult;
}

bool RAbundVectorTestFixture::TestSet(const int binNumber, const int newBinSize, const bool expectedResult) {
    Setup();
    const int currentRankSize = radbund.getMaxRank();
    radbund.set(binNumber, newBinSize);
    const bool result = radbund.getMaxRank() != currentRankSize;
    TearDown();
    return expectedResult == result;

}

bool RAbundVectorTestFixture::TestGet(const int index, const int expectedResult) {
    Setup();
    const bool result = radbund.get(index) == expectedResult;
    TearDown();
    return result;
}

bool RAbundVectorTestFixture::TestGetOverload(const std::vector<int>& expectedResult) {
    Setup();
    const auto result = radbund.get();
    TearDown();
    return result.size() == expectedResult.size();
}

bool RAbundVectorTestFixture::TestPushBack(const int binSize, const bool expectedResult) {
    Setup();
    const int currentBinSize = radbund.getNumBins();
    radbund.push_back(binSize);
    const auto result = radbund.getNumBins() != currentBinSize;
    TearDown();
    return result == expectedResult;
}

bool RAbundVectorTestFixture::TestResize(const int size, const bool expectedResult) {
    Setup();
    radbund.resize(size);
    const bool result = radbund.size() == size;
    TearDown();
    return result == expectedResult;
}


bool RAbundVectorTestFixture::TestSize(const int expectedResult) {
    Setup();
    const bool result = radbund.size() == expectedResult;
    TearDown();
    return result;
}

bool RAbundVectorTestFixture::TestQuicksort(const bool expectedResult) {
    Setup();
    const bool result = radbund.quicksort();
    TearDown();
    return result == expectedResult;
}

bool RAbundVectorTestFixture::TestClear(const bool expectedResult) {
    Setup();
    radbund.clear();
    const bool result = radbund.size() <= 0;
    TearDown();
    return result == expectedResult;
}

bool RAbundVectorTestFixture::TestGetNumBins(const int expectedResult) {
    Setup();
    const bool result = radbund.getNumBins() == expectedResult;
    TearDown();
    return result;
}

bool RAbundVectorTestFixture::TestGetNumSeqs(const int expectedResult) {
    Setup();
    const bool result = radbund.getNumSeqs() == expectedResult;
    TearDown();
    return result;
}

bool RAbundVectorTestFixture::TestGetMaxRank(const int expectedResult) {
    Setup();
    const bool result = radbund.getMaxRank() == expectedResult;
    TearDown();
    return result;
}

void RAbundVectorTestFixture::Setup() {
    const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
    const std::vector<double> total{10, 20, 30, 40, 50, 60};
    const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
        Rcpp::Named("Representative Sequence") = compounds,
        Rcpp::Named("total") = total,
        Rcpp::Named("nogroup") = total);
    CountTableAdapter countTable;
    countTable.CreateDataFrameMap(dataframe);
    MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25},
        0.2, false, countTable);
    // adapter.CreateSparseMatrix();
    radbund = adapter.CreateListVector().getRAbundVector();
}

void RAbundVectorTestFixture::TearDown() {
    radbund = RAbundVector();
}
