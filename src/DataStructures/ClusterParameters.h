//
// Created by Gregory Johnson on 9/9/26.
//

#ifndef REFACTOR_CLUSTERPARAMETERS_H
#define REFACTOR_CLUSTERPARAMETERS_H
#include <string>
#include <unordered_map>
#include <utility>

class ClusterParameters {
public:
    ~ClusterParameters() = default;
    ClusterParameters() = default;
    explicit ClusterParameters(std::string  method):clusterMethod(std::move(method)) {}
     std::string GetClusterType() {return clusterMethod; }
private:
    std::string clusterMethod;
    std::unordered_map<std::string, std::string> parameters;
};
#endif //REFACTOR_CLUSTERPARAMETERS_H