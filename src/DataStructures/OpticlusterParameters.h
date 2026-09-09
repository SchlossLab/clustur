//
// Created by Gregory Johnson on 9/9/26.
//

#ifndef REFACTOR_OPTICLUSTERPARAMETERS_H
#define REFACTOR_OPTICLUSTERPARAMETERS_H
#include "ClusterParameters.h"

class OptiClusterParameters final : ClusterParameters {
public:
    std::string GetClusterType() override;
    OptiClusterParameters(const OptiClusterParameters&) = default;
    OptiClusterParameters(OptiClusterParameters&&) = default;

};
#endif //REFACTOR_OPTICLUSTERPARAMETERS_H