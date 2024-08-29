//
// Created by Gregory Johnson on 8/21/24.
//

#include "MothurDependencies/OpticlusterData.h"

void OpticlusterData::AddToData(OptiClusterInformation& information) {
    opticlusterInformation.emplace_back(information);
}
