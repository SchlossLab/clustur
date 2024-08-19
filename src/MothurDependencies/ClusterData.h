//
// Created by Gregory Johnson on 8/19/24.
//

#ifndef CLUSTERDATA_H
#define CLUSTERDATA_H
#include <string>
#include <unordered_map>
#include <vector>

struct ClusterData {
    std::string label;
    int numberOfOtu;
    std::unordered_map<int, std::vector<int>> clusters;
};
#endif //CLUSTERDATA_H
