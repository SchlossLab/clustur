//
// Created by Gregory Johnson on 8/19/24.
//

#ifndef CLUSTERDATA_H
#define CLUSTERDATA_H
#include <string>
#include <unordered_map>
#include <vector>

#include "ClusterExport.h"

struct ClusterInformation {
    std::string label;
    int numberOfOtu;
    std::string clusterBins;
};

class ClusterData final : public ClusterExport {
public:
    explicit ClusterData(std::string header):headers(std::move(header)){}
    std::string Print() const override;
    void AddToData(ClusterInformation&);
    std::string GetLabel(int) const override;

private:
    std::vector<ClusterInformation> clusterInformation;
    std::string headers;
};

#endif //CLUSTERDATA_H
