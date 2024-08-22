//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef SHAREDFILE_H
#define SHAREDFILE_H
#include <unordered_map>
#include <vector>
#include <fstream>
#include "OtuAbundancePair.h"




class SharedFile {
public:
    explicit SharedFile(const std::unordered_map<std::string, std::vector<OTUAbundancePair>>& sharedMap,
        const int binSize):
    otuAbundance(sharedMap),
    largestBin(binSize){}
    void PrintData() const;

private:
    std::unordered_map<std::string, std::vector<OTUAbundancePair>> otuAbundance;
    int largestBin;
};



#endif //SHAREDFILE_H
