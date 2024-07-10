//
// Created by gregj on 7/9/2024.
//

#ifndef CLUSTERTESTFIXTURE_H
#define CLUSTERTESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/Cluster.h"


class ClusterTestFixture : public TestFixture {
public:
    ClusterTestFixture() {};
    bool TestUpdateClusterWorks();
    bool TestGetTagGetter();
    bool TestSetWantedMap();
    bool TestGetSeqToBin();
    bool TestUpdateDistanceOverload();
    bool TestClusterBins();
    bool TestClusterNames();
    bool TestUpdateMap();
    ~ClusterTestFixture() override;

private:
    void Setup() override;
    void TearDown() override;
    Cluster* cluster;
};



#endif //CLUSTERTESTFIXTURE_H
