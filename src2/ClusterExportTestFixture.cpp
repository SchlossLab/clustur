//
// Created by Gregory Johnson on 8/29/24.
//

#include "Tests/ClusterExportTestFixture.h"

bool ClusterExportTestFixture::TestSetListVector(const ListVector& listVector, const std::string& label, const bool result) {
    Setup();
    clusterExport->SetListVector(listVector, label);
    const auto pair = clusterExport->GetListVector();
    const std::string otherLabel = pair.label;
    const int otherSize = pair.listVector.size();
    const bool finalResult = result == (otherLabel == label && otherSize == listVector.size());
    TearDown();
    return finalResult;
}

bool ClusterExportTestFixture::TestGetListVector(const bool result) {
    Setup();
    const ListVector vector;
    clusterExport->SetListVector(vector, "0.00");
    const auto pair = clusterExport->GetListVector();
    const int otherSize = pair.listVector.size();
    const bool finalResult = result == (otherSize == vector.size());
    TearDown();
    return finalResult;
}

void ClusterExportTestFixture::Setup() {
    clusterExport = new ClusterExport();
}

void ClusterExportTestFixture::TearDown() {
    delete clusterExport;
}
