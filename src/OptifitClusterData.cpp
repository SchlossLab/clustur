//
// Created by Gregory Johnson on 7/24/26.
//

#include "DataExporters/OptifitClusterData.h"

bool OptifitClusterData::AddToData(OptifitClusterInformation& information) {
    if(information.label.empty())
        return false;
    optiClusterInformation.emplace_back(information);
    return true;
}
