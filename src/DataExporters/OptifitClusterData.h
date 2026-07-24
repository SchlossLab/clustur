//
// Created by Gregory Johnson on 7/24/26.
//

#ifndef REFACTOR_OPTIFITCLUSTERDATA_H
#define REFACTOR_OPTIFITCLUSTERDATA_H
#include "ClusterExport.h"


struct OptifitClusterInformation {
    int numberOfOtu;
    std::string clusterBins;
    std::string label;
};

class OptifitClusterData final : public ClusterExport {
public:
    explicit OptifitClusterData(std::string header):headers(std::move(header)){}
    bool AddToData(OptifitClusterInformation&);

private:
    std::vector<OptifitClusterInformation> optiClusterInformation;
    std::string headers;
};


#endif //REFACTOR_OPTIFITCLUSTERDATA_H