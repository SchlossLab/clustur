//
// Created by gregj on 8/20/2024.
//

#include "MothurDependencies/ClusterData.h"

#include "TestHelpers/TestHelper.h"

std::string ClusterData::Print() const {
    std::string results = headers;
    for(const auto& cluster : clusterInformation) {
        results += cluster.label + "\t" + cluster.clusterBins;
    }
    return results;
}

void ClusterData::AddToData(ClusterInformation& information) {
    clusterInformation.emplace_back(information);
}
