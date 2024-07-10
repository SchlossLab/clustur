//
// Created by gregj on 7/9/2024.
//

#include "ClusterTestFixture.h"

bool ClusterTestFixture::TestGetTagGetter(Cluster* clust, const std::string& expectedResult) {
    Setup();
    cluster = clust;
    const std::string tag = cluster->getTag();
    TearDown();
    return tag == expectedResult;
}

bool ClusterTestFixture::TestGetSeqToBin(Cluster* clust, const std::map<std::string, int>& expectedMap) {
    Setup();
    cluster = clust;
    const auto bin = cluster->getSeqtoBin();
    TearDown();
    return bin == expectedMap;
}

bool ClusterTestFixture::TestUpdateDistanceOverload(Cluster* clust, PDistCell& row, PDistCell& col,
    const bool expected) {
    Setup();
    cluster = clust;
    const bool result = clust->updateDistance(row, col);
    TearDown();
    return result == expected;
}

bool ClusterTestFixture::TestClusterBins(RAbundVector* rAbdund, const bool expectedResult) {
    Setup();
    rabund = rAbdund;
    const bool result = clusterBins();
    TearDown();
    return result == expectedResult;
}

bool ClusterTestFixture::TestClusterNames(ListVector* listVector, const bool expectedResult) {
    Setup();
    list = listVector;
    const bool result = clusterNames();
    TearDown();
    return result == expectedResult;
}

bool ClusterTestFixture::TestUpdateMap(ListVector* listVector, const bool expectedResult) {

    Setup();
    list = listVector;
    const bool result = updateMap();
    TearDown();
    return result == expectedResult;
}

ClusterTestFixture::~ClusterTestFixture() {
    delete cluster;
}

void ClusterTestFixture::Setup() {
    cluster = nullptr;
}

void ClusterTestFixture::TearDown() {
    delete cluster;
}
