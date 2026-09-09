//
// Created by Gregory Johnson on 9/9/26.
//

#ifndef REFACTOR_CLUSTERPARAMETERS_H
#define REFACTOR_CLUSTERPARAMETERS_H
#include <string>
class ClusterParameters {
public:
    virtual ~ClusterParameters() = default;
    ClusterParameters() = default;
    virtual std::string GetClusterType();
};
#endif //REFACTOR_CLUSTERPARAMETERS_H