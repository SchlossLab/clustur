//
// Created by Gregory Johnson on 8/19/24.
//

#ifndef CLUSTERRESULT_H
#define CLUSTERRESULT_H
#include "ClusterData.h"
#include <string>
#include <vector>
#include <unordered_map>

class ClusterResult {
public:
    std::string ResultsToString() const;
    std::vector<ClusterData> clusterData;
};

#endif //CLUSTERRESULT_H
