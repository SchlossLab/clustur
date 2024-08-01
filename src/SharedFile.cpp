//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFile.h"

void SharedFile::PrintData() const {
    std::string headers;
    std::string data;
    for(const auto& abundancePairs : otuAbundance) {
        headers += abundancePairs.otu + '\t';
        data += std::to_string(abundancePairs.abundance) + '\t';
    }
    data[data.size() - 1] = ' ';
    std::ofstream sharedFileOutput("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/SharedFile.txt");
    sharedFileOutput << headers + "\n" + data + "\n";
    sharedFileOutput.close();
}
