//
// Created by gregj on 8/20/2024.
//

#include "MothurDependencies/ClusterData.h"

#include "TestHelpers/TestHelper.h"

void ClusterData::AddToData(ClusterInformation& information) {
    clusterInformation.emplace_back(information);
}