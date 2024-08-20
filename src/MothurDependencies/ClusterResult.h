//
// Created by Gregory Johnson on 8/19/24.
//

#ifndef CLUSTERRESULT_H
#define CLUSTERRESULT_H
#include "ClusterData.h"
#include <string>
#include <utility>
#include <vector>

class ClusterResult {
public:
    ClusterResult() = default;
    std::string ResultsToString(const ClusterExport&) const;
};

#endif //CLUSTERRESULT_H
