//
// Created by Gregory Johnson on 8/29/24.
//

#include "Tests/ClusterExportTestFixture.h"

bool ClusterExportTestFixture::TestSetListVector(ListVector& listVector, const std::string& label, const bool result) {
    Setup();
    clusteExport->SetListVector(listVector, label);
    const auto pair = clusteExport->GetListVector();
    const std::string otherLabel = pair.label;
    const int otherSize = pair.listVector->size();
    const bool finalResult = result == (otherLabel == label && otherSize == listVector.size());
    TearDown();
    return finalResult;
}

bool ClusterExportTestFixture::TestGetListVector(const bool result) {
    Setup();
    auto* vector = new ListVector();
    clusteExport->SetListVector(*vector, "0.00");
    const auto pair = clusteExport->GetListVector();
    const int otherSize = pair.listVector->size();
    const bool finalResult = result == (otherSize == vector->size());
    TearDown();
    return finalResult;
}

void ClusterExportTestFixture::Setup() {
    clusteExport = new ClusterExport();
}

void ClusterExportTestFixture::TearDown() {
    delete(clusteExport);
}
