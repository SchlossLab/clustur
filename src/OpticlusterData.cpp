//
// Created by Gregory Johnson on 8/21/24.
//

#include "MothurDependencies/OpticlusterData.h"

void OpticlusterData::AddToData(OptiClusterInformation& information) {
    opticlusterInformation.emplace_back(information);
}

std::string OpticlusterData::Print() const {
    std::string results = headers;
    for(const auto& cluster : opticlusterInformation) {
        results += cluster.clusterBins;
    }
    return results;


}

std::string OpticlusterData::GetLabel(const int index) const {
    return opticlusterInformation[index].label;
}
