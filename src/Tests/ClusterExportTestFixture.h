//
// Created by Gregory Johnson on 8/29/24.
//

#ifndef CLUSTEREXPORTTESTFIXTURE_H
#define CLUSTEREXPORTTESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/ClusterExport.h"
#include "../MothurDependencies/ListVector.h"


class ClusterExportTestFixture final : public TestFixture {
public:
    bool TestSetListVector(const ListVector &listVector, const std::string &label, bool result);
    bool TestGetListVector(bool result);
private:
    void Setup() override;
    void TearDown() override;
    ClusterExport* clusterExport = nullptr;
};



#endif //CLUSTEREXPORTTESTFIXTURE_H
