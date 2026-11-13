//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_OPTICLUSTERDATA_H
#define REFACTOR_OPTICLUSTERDATA_H

#include "ClusterExport.h"
#include <string>
#include <vector>


struct OptiClusterInformation {
    int numberOfOtu;
    std::string clusterBins;
    std::string label;
};

class OptiClusterData final: public ClusterExport {
public:
    explicit OptiClusterData(std::string header):headers(std::move(header)){}
    bool AddToData(OptiClusterInformation&);

private:
    std::vector<OptiClusterInformation> optiClusterInformation;
    std::string headers;
};
#endif //REFACTOR_OPTICLUSTERDATA_H