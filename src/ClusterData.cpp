//
// Created by gregj on 8/20/2024.
//

#include "MothurDependencies/ClusterData.h"

#include "TestHelpers/TestHelper.h"

std::string ClusterData::Print() const {
    std::string results = "label\t";
   
    const int largestBinSize = GetLargestBinSize();
    const int length = std::to_string(largestBinSize).length();
    for(int i = 0; i < GetLargestBinSize(); i++) {
        std::stringstream ss;
        ss << std::setw(length) << std::setfill('0') << (i + 1);
        results += "OTU" + ss.str() + "\t";
    }
    results[results.length() - 1] = ' ';
    results += "\n";
    for(const auto& cluster : clusterInformation) {
        results += cluster.label + cluster.clusterBins;
    }
    return results;
}

void ClusterData::AddToData(ClusterInformation& information) {
    clusterInformation.emplace_back(information);
}

std::string ClusterData::GetLabel(const int index) const {
    return clusterInformation[index].label;
}
