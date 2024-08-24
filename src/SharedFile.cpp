//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFile.h"

#include "MothurDependencies/ClusterExport.h"

// TODO: Fix this function so that it exports a clean shared file
void SharedFile::PrintData() const {
    std::string headers = "label\t";
    for(int i = 0; i < largestBin; i++) {
        headers += "OTU" + std::to_string(i) + "\t";
    }
    std::string sharedFile;
    for(const auto& abund : otuAbundance) {
        // const size_t size = abund.second.size(); bin size
        std::string data = abund.first + "\t"; // label
        for(const auto& otuPairs : abund.second) {
            data += std::to_string(otuPairs.abundance) + '\t';
        }
        sharedFile += data + "\n";
    }
    sharedFile[sharedFile.size() - 1] = ' ';
    std::ofstream sharedFileOutput("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/SharedFile.txt");
    sharedFileOutput << headers + "\n" + sharedFile + "\n";
    sharedFileOutput.close();
}
