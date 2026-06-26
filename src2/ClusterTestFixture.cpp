//
// Created by gregj on 7/9/2024.
//

#include "Tests/ClusterTestFixture.h"

bool ClusterTestFixture::TestGetTagGetter(Cluster* clust, const std::string& expectedResult) {
    Setup();
    const std::string tag = clust->getTag();
    TearDown();
    return tag == expectedResult;
}

bool ClusterTestFixture::TestGetSeqToBin(Cluster* clust, const std::map<std::string, int>& expectedMap) {
    Setup();
    const auto bin = clust->getSeqtoBin();
    TearDown();
    return bin == expectedMap;
}

bool ClusterTestFixture::TestClusterBins(RAbundVector* rAbdund, SparseDistanceMatrix* distanceMatrix,
    const bool expectedResult) {
    Setup();
    dMatrix = distanceMatrix;
    smallCol = dMatrix->getSmallestCell(smallRow);
    rabund = rAbdund;
    const int currentMaxRank = rAbdund->getMaxRank();
    clusterBins();
    const bool result = rAbdund->getMaxRank() != currentMaxRank;
    TearDown();
    return result == expectedResult;
}

bool ClusterTestFixture::TestClusterNames(ListVector* listVector, SparseDistanceMatrix* distanceMatrix,
    const bool expectedResult) {
    Setup();
    dMatrix = distanceMatrix;
    list = listVector;
    smallCol = dMatrix->getSmallestCell(smallRow);
    const bool result = clusterNames();
    TearDown();
    return result == expectedResult;
}

bool ClusterTestFixture::TestUpdateMap(ListVector* listVector, SparseDistanceMatrix* distanceMatrix, const bool expectedResult) {

    Setup();
    list = listVector;
    dMatrix = distanceMatrix;
    smallCol = dMatrix->getSmallestCell(smallRow);
    const bool result = updateMap();
    TearDown();
    return result == expectedResult;
}


void ClusterTestFixture::Setup() {
    cluster = nullptr;
    rabund = nullptr;
    dMatrix = nullptr;
}

void ClusterTestFixture::TearDown() {
    delete cluster;
}
