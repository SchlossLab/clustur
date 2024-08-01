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
    explicit SharedFile(const std::vector<OTUAbundancePair>& otuAbundance):otuAbundance(otuAbundance) {}
    void PrintData() const;

private:
    std::vector<OTUAbundancePair> otuAbundance;
};



#endif //SHAREDFILE_H
