//
// Created by Gregory Johnson on 6/9/26.
//

#include "OptiClusterData.h"

bool OptiClusterData::AddToData(OptiClusterInformation& information) {
    if(information.label.empty())
        return false;
    optiClusterInformation.emplace_back(information);
    return true;
}
