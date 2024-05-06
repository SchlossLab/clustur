//
// Created by Gregory Johnson on 5/3/24.
//

#include "tests/OptiClusterTestFixture.h"

bool OptiClusterTestFixture::OptiClusterInitalizesProperly(double metricValue, bool randomize,
    const std::string &initializationString, int expectedResult) {
}

bool OptiClusterTestFixture::OptiClusterUpdatesCorrectly(double metric, bool expectedResult) {
}

bool OptiClusterTestFixture::OpticlusterReturnsTheCorrectStats(double tp, double tn, double fp, double fn,
    std::vector<double> expectedResults) {
}

bool OptiClusterTestFixture::OptiClusterGetsTheCorrectCloseFarCounts(long long sequence, long long newBin,
    std::vector<double> expectedResults) {
}

bool OptiClusterTestFixture::OptiClusterReturnsTheCorrectList(ListVector *expectedResult) {
}

bool OptiClusterTestFixture::OptiClusterGetsTheCorrectTag(const std::string &expectedResult) {
}

bool OptiClusterTestFixture::OptiClusterGetsTheCorrectNumberOfBins(long long expectedResult) {
}

bool OptiClusterTestFixture::OptiClusterFindsCorrectInsert(long long expectedResult) {
}

void OptiClusterTestFixture::Setup() {
}

void OptiClusterTestFixture::TearDown() {
}
