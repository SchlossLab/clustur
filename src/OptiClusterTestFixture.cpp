//
// Created by Gregory Johnson on 5/3/24.
//

#include "tests/OptiClusterTestFixture.h"

#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/Metrics/mcc.hpp"


bool OptiClusterTestFixture::OptiClusterInitalizesProperly(double metricValue, const bool randomize,
                                                           const std::string &initializationString,
                                                           const int expectedResult) {
    Setup();
    const int returnedValue = optiCluster->initialize(metricValue, randomize,initializationString);
    TearDown();
    return returnedValue == expectedResult;
}

bool OptiClusterTestFixture::OptiClusterUpdatesCorrectly(double metric, const bool expectedResult) {
    Setup();
    const bool result = optiCluster->update(metric);
    TearDown();
    return result == expectedResult;
}

bool OptiClusterTestFixture::OpticlusterReturnsTheCorrectStats(double tp, double tn, double fp, double fn,
    const std::vector<double>& expectedResults) {
    Setup();
    const auto stats = optiCluster->getStats(tp, tn, fp, fn);
    TearDown();
    return stats == expectedResults;
}

bool OptiClusterTestFixture::OptiClusterGetsTheCorrectCloseFarCounts(const long long sequence, const long long newBin,
    const std::vector<double>& expectedResults) {
    Setup();
    double metric = 0;
    optiCluster->initialize(metric, false, "");
    const auto farCounts = optiCluster->getCloseFarCounts(sequence, newBin);
    TearDown();
    const bool val = farCounts == expectedResults;
    return val;
}

bool OptiClusterTestFixture::OptiClusterReturnsTheCorrectList(const double expectedDataSize) {
    Setup();
    const auto list = optiCluster->getList();
    TearDown();
    return list->getNumBins() == expectedDataSize;
}

bool OptiClusterTestFixture::OptiClusterGetsTheCorrectTag(const std::string &expectedResult) {
    Setup();
    const auto tags = optiCluster->getTag();
    TearDown();
    return tags == expectedResult;
}

bool OptiClusterTestFixture::OptiClusterGetsTheCorrectNumberOfBins(const long long expectedResult) {
    Setup();
    const long long bins = optiCluster->getNumBins();
    TearDown();
    return bins == expectedResult;
}


void OptiClusterTestFixture::Setup() {
    const auto xVals = std::vector<int>{0,0,0,1,1,2,3};
    const auto yVals = std::vector<int>{1,2,4,2,4,4,4};
    const auto data = std::vector<double>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
    OptimatrixAdapter adapter(0.03);
    OptiData *matrix = adapter.ConvertToOptimatrix(xVals, yVals, data);
    optiCluster = std::make_unique<OptiCluster>(matrix, new MCC(), 0);;
}

void OptiClusterTestFixture::TearDown() {
    optiCluster.reset(nullptr);
}
