//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFileBuilder.h"

#include "MothurDependencies/ClusterExport.h"

// TODO Comment this code
SharedFile* SharedFileBuilder::BuildSharedFile(const std::unordered_map<std::string, ListVector> &listVectorMap,
    const CountTableAdapter &countTable) {
    Utils utils;
    double largestCutoff = 0;
    std::string largestCutoffLabel;
    for(auto& map : listVectorMap) {
        const double cutoffLabel = std::stod(map.first);
        if(largestCutoff < cutoffLabel) {
            largestCutoff = cutoffLabel;
            largestCutoffLabel = map.first;
        }
    }
    std::vector<SharedAbundance> abundancesList;
    std::vector<std::string> groups = countTable.GetGroups();
    ListVector correctListVector = listVectorMap.at(largestCutoffLabel);
    const int seqs = correctListVector.getNumSeqs();
    int count = 1;
    for(int i = 0; i < seqs; i++) {
        std::string samples = correctListVector.get(i);
        if(samples.empty())
            continue;
        std::vector<std::string> splitSamples;
        std::string otuName = "OTU" + std::to_string(count++);
        utils.splitAtComma(samples, splitSamples);
        std::unordered_map<std::string, double> totalAbundanceInEachGroup;
        for(const auto& sample : splitSamples) {
            for(const auto& group : groups) { // Its already in tidy form
                totalAbundanceInEachGroup[group] += countTable.FindAbundanceBasedOnGroup(group, sample);
            }
        }
        for(const auto& groupTotals : totalAbundanceInEachGroup) {
            abundancesList.emplace_back(groupTotals.first, otuName,
                 largestCutoffLabel, groupTotals.second);
        }
    }
    // Have to make it group totals, not the abundance of each sample at each group
   return new SharedFile(abundancesList);

}
