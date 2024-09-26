//
// Created by gregj on 8/20/2024.
//

#include "MothurDependencies/ClusterData.h"

#include "TestHelpers/TestHelper.h"

bool ClusterData::AddToData(ClusterInformation& information) {
    if(information.label.empty())
        return false;
    clusterInformation.emplace_back(information);
    return true;
}
