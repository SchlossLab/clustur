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
    explicit ClusterParameters(std::string method):clusterMethod(std::move(method)) {
        parameters["iters"] = "100";
        parameters["delta"] = "1";
        parameters["precision"] = "100";
        parameters["initialize"] = "singleton";
    }
    [[nodiscard]] const std::string& GetClusterType() const {return clusterMethod; }
    [[nodiscard]] const std::unordered_map<std::string, std::string>& GetClusterParameters() const {return parameters;}
    void SetParameters(const std::string& key, const std::string& value);
private:
    std::string clusterMethod;
    std::unordered_map<std::string, std::string> parameters;
};
#endif //REFACTOR_CLUSTERPARAMETERS_H