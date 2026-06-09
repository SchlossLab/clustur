//
// Created by Gregory Johnson on 6/9/26.
//

#include "ClusterData.h"

bool ClusterData::AddToData(ClusterInformation& information) {
    if(information.label.empty())
        return false;
    clusterInformation.emplace_back(information);
    return true;
}
