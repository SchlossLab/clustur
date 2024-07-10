//
// Created by gregj on 7/9/2024.
//

#include "ClusterTestFixture.h"

bool ClusterTestFixture::TestUpdateClusterWorks() {
}

bool ClusterTestFixture::TestGetTagGetter() {
}

bool ClusterTestFixture::TestSetWantedMap() {
}

bool ClusterTestFixture::TestGetSeqToBin() {
}

bool ClusterTestFixture::TestUpdateDistanceOverload() {
}

bool ClusterTestFixture::TestClusterBins() {
}

bool ClusterTestFixture::TestClusterNames() {
}

bool ClusterTestFixture::TestUpdateMap() {
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
