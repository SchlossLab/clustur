//
// Created by Gregory Johnson on 8/21/24.
//

#include "MothurDependencies/OpticlusterData.h"

bool OpticlusterData::AddToData(OptiClusterInformation& information) {
    if(information.label.empty())
        return false;
    opticlusterInformation.emplace_back(information);
    return true;
}
