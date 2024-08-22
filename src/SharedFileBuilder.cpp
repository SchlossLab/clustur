//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFileBuilder.h"

#include "MothurDependencies/ClusterExport.h"


//TODO: Find a way to incorporate groups
//TODO Find a way to incorporate the map correctly
SharedFile* SharedFileBuilder::BuildSharedFile(const std::unordered_map<std::string, RAbundVector>& rAdbundMap) {
    std::unordered_map<std::string, std::vector<OTUAbundancePair>> sharedFileMap;
    int largestBin = 0;
    int count = 0;
    for(const auto& abundVectorMap : rAdbundMap) {
        std::string label = abundVectorMap.first;
        auto rAbound = abundVectorMap.second;
        const int vecSize = rAbound.getNumBins();
        for(int i = 0; i < vecSize; i++) 
        {
            const int abundance = rAbound.get(i);
            if(abundance <= 0)
                continue;
            sharedFileMap[label].emplace_back(abundance,
           "OTU" + std::to_string(i + 1),"nogroup", label);
        }
         if(vecSize >  sharedFileMap[label].size()) {
            largestBin = sharedFileMap[label].size();
        }
    }
    return new SharedFile(sharedFileMap, largestBin);


    //
    //
    //
    // std::vector<OTUAbundancePair> abundanceMap(size);
    // std::vector<OTUAbundancePair> pair;
    // for(int i = 0; i < size; i++) {
    //     std::vector<std::string> container;
    //     const int abundance = rAdbundVector->get(i);
    //     abundanceMap[i] = OTUAbundancePair(abundance,
    //         "OTU" + std::to_string(i + 1),"nogroup", clusterExport->GetLabel(i));
    //     // Get abundance at each otu level
    // }
    // return new SharedFile(abundanceMap);
}
