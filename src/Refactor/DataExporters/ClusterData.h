//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_CLUSTERDATA_H
#define REFACTOR_CLUSTERDATA_H

#include <string>
#include <unordered_map>
#include <vector>

#include "ClusterExport.h"

struct ClusterInformation {
    std::string label;
    int numberOfOtu;
    std::string clusterBins;
};

class ClusterData final : public ClusterExport {
public:
    explicit ClusterData(std::string header):headers(std::move(header)){}
    bool AddToData(ClusterInformation&);

private:
    std::vector<ClusterInformation> clusterInformation;
    std::string headers;
};


#endif //REFACTOR_CLUSTERDATA_H