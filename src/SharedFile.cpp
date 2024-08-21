//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFile.h"

#include "MothurDependencies/ClusterExport.h"

// TODO: Fix this function so that it exports a clean shared file
void SharedFile::PrintData(ClusterExport *clusterExport) const {
    std::string headers = "label\t";
    std::string sharedFile;
    for(size_t i = 0; i < clusterExport->GetAbundanceCount(); i++) {
        std::string data = otuAbundance[i].label + "\t";
        for (const auto &abundancePairs: otuAbundance) {
            //headers += abundancePairs.otu + '\t';
            data += std::to_string(abundancePairs.abundance) + '\t';
        }
    }
    sharedFile[sharedFile.size() - 1] = ' ';
    std::ofstream sharedFileOutput("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/SharedFile.txt");
    sharedFileOutput << headers + "\n" + sharedFile + "\n";
    sharedFileOutput.close();
}
