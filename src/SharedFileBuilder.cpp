//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFileBuilder.h"

#include "MothurDependencies/ClusterExport.h"


//TODO: Find a way to incorporate groups
SharedFile* SharedFileBuilder::BuildSharedFile(RAbundVector *rAdbundVector, const ClusterExport* clusterExport) {
    const int size = rAdbundVector->size();


    std::vector<OTUAbundancePair> abundanceMap(size);
    std::vector<OTUAbundancePair> pair;
    for(int i = 0; i < size; i++) {
        std::vector<std::string> container;
        const int abundance = rAdbundVector->get(i);
        abundanceMap[i] = OTUAbundancePair(abundance,
            "OTU" + std::to_string(i + 1),"nogroup", clusterExport->GetLabel(i));
        // Get abundance at each otu level
    }
    return new SharedFile(abundanceMap);
}
