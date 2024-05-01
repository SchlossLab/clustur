//
// Created by Gregory Johnson on 5/1/24.
//

#include "tests/ClusterCommandTestFixture.h"

bool ClusterCommandTestFixture::TestRunOptiClusterReturnsCorrectData(const OptiMatrix *mockMatrix,
    const std::string& expectedResult) {

}

bool ClusterCommandTestFixture::TestSetIterationsWorksCorrectly(const int iterations, const int expectedIterations) {

}

void ClusterCommandTestFixture::Setup() {
    clusterCommand = std::make_unique<ClusterCommand>();
}

void ClusterCommandTestFixture::TearDown() {
    clusterCommand.reset();
}
