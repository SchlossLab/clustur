//
// Created by Gregory Johnson on 5/3/24.
//

#ifndef OPTICLUSTERTESTFIXTURE_H
#define OPTICLUSTERTESTFIXTURE_H
#include <memory>
#include "TestFixture.h"
#include "../MothurDependencies/OptiCluster.h"


class OptiClusterTestFixture : TestFixture {
    bool OptiClusterInitalizesProperly(double metricValue, bool randomize,const std::string& initializationString,
        int expectedResult);
    bool OptiClusterUpdatesCorrectly(double metric, bool expectedResult);
    bool OpticlusterReturnsTheCorrectStats(double tp, double tn, double fp, double fn,
        std::vector<double> expectedResults);
    bool OptiClusterGetsTheCorrectCloseFarCounts(long long sequence, long long newBin,
        std::vector<double> expectedResults);
    bool OptiClusterReturnsTheCorrectList(ListVector* expectedResult);
    bool OptiClusterGetsTheCorrectTag(const std::string& expectedResult);
    bool OptiClusterGetsTheCorrectNumberOfBins(long long expectedResult);
    bool OptiClusterFindsCorrectInsert(long long expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    std::unique_ptr<OptiCluster> optiCluster;
};



#endif //OPTICLUSTERTESTFIXTURE_H
