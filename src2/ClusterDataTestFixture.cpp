//
// Created by Gregory Johnson on 8/29/24.
//

#include "Tests/ClusterDataTestFixture.h"

bool ClusterDataTestFixture::CheckAddToData(ClusterInformation& information, const bool expectedResult) {
    Setup();
    const bool result = data->AddToData(information);
    TearDown();
    return expectedResult == result;
}

void ClusterDataTestFixture::Setup() {
    data = new ClusterData("");
}

void ClusterDataTestFixture::TearDown() {
    delete(data);
}
