//
// Created by gregj on 7/9/2024.
//

#ifndef CLUSTERTESTFIXTURE_H
#define CLUSTERTESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/Cluster.h"


class ClusterTestFixture : public TestFixture , public Cluster {
public:
    ClusterTestFixture() = default;
    bool TestGetTagGetter(Cluster *clust, const std::string &expectedResult);
    bool TestGetSeqToBin(Cluster* clust, const std::map<std::string, int>& expectedMap);
    bool TestUpdateDistanceOverload(Cluster* clust, PDistCell& row, PDistCell& col, bool expected);
    bool TestClusterBins(RAbundVector* rAbdund, bool expectedResult);
    bool TestClusterNames(ListVector* listVector, bool expectedResult);
    bool TestUpdateMap(ListVector* listVector, bool expectedResult);
    ~ClusterTestFixture() override;

private:
    void Setup() override;
    void TearDown() override;
    Cluster* cluster = nullptr;
};



#endif //CLUSTERTESTFIXTURE_H
