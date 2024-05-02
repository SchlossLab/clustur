//
// Created by Gregory Johnson on 5/1/24.
//

#include "tests/ClusterCommandTestFixture.h"

bool ClusterCommandTestFixture::TestRunOptiClusterReturnsCorrectData(OptiMatrix *mockMatrix,
    const std::string& expectedResult) {
    Setup();
    const std::string result = clusterCommand->runOptiCluster(mockMatrix);
    TearDown();
    return result == expectedResult;
}

bool ClusterCommandTestFixture::TestSetIterationsWorksCorrectly(const int iterations, const bool expectResult) {
    Setup();
    const bool result = clusterCommand->SetMaxIterations(iterations);
    TearDown();
    return result == expectResult;
}

void ClusterCommandTestFixture::Setup() {
    clusterCommand = std::make_unique<ClusterCommand>();
}

void ClusterCommandTestFixture::TearDown() {
    clusterCommand.reset();
}
