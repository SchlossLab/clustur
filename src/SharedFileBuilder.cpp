//
// Created by Gregory Johnson on 7/31/24.
//

#include "SharedFileData/SharedFileBuilder.h"
#include "DataExporters/ClusterExport.h"

// TODO Comment this code
// TODO We may need to build a traditional file builder...So we can output a dataframe of how the clusters are (list)
SharedFile* SharedFileBuilder::BuildSharedFile(const ListVector &listVector,
                                               const CountTableAdapter &countTable, const std::string &binName) {
    std::string largestCutoffLabel = listVector.getLabel();
    std::vector<SharedAbundance> abundancesList;
    const std::vector<std::string> groups = countTable.GetGroups();
    const int seqs = listVector.getNumSeqs();
    int count = 1;
    for(int i = 0; i < seqs; i++) { // O(N^3) -> Lets make this O(N^2) or O(Nlog(N))
        std::string samples = listVector.get(i);
        if(samples.empty())
            continue;
        std::vector<std::string> splitSamples;
        std::string otuName = binName + std::to_string(count++);
        Utils::splitAtComma(samples, splitSamples);
        std::unordered_map<std::string, double> totalAbundanceInEachGroup;
        for(const auto& sample : splitSamples) {
            for(const auto& group : groups) { // Its already in tidy form
                totalAbundanceInEachGroup[group] += countTable.FindAbundanceBasedOnGroup(group, sample);
            }
        }
        for(const auto&[fst, snd] : totalAbundanceInEachGroup) {
            abundancesList.emplace_back(fst, otuName,
                 largestCutoffLabel, snd);
        }
    }
    
   return new SharedFile(abundancesList);
}
