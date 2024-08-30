//
// Created by Gregory Johnson on 8/29/24.
//

#include "Tests/OptiClusterDataTestFixture.h"

bool OptiClusterDataTestFixture::CheckAddToData(OptiClusterInformation& information, bool expectedResult) {
    Setup();
    const bool result = data->AddToData(information);
    TearDown();
    return expectedResult == result;
}

void OptiClusterDataTestFixture::Setup() {
    data = new OpticlusterData("", 0);
}

void OptiClusterDataTestFixture::TearDown() {
    delete(data);
}
