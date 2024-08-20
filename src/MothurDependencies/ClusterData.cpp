//
// Created by gregj on 8/20/2024.
//

#include "ClusterData.h"

std::string ClusterData::Print() const {
    std::string results = headers;
    for(const auto& cluster : clusterInformation) {
        results += cluster.label + std::to_string(cluster.numberOfOtu) + cluster.clusterBins + "\n";
    }
    return results;
}

void ClusterData::AddToData(ClusterInformation& information) {
    clusterInformation.emplace_back(information);
}
