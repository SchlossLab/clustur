//
// Created by Gregory Johnson on 8/21/24.
//

#ifndef OPTICLUSTERDATA_H
#define OPTICLUSTERDATA_H
#include "ClusterExport.h"
#include <string>
#include <vector>

struct OptiClusterInformation {
    int numberOfOtu;
    std::string clusterBins;
    std::string label;
};

class OpticlusterData final: public ClusterExport{
public:
    explicit OpticlusterData(std::string header, const double cutoff):headers(std::move(header)){}
    bool AddToData(OptiClusterInformation&);

private:
    std::vector<OptiClusterInformation> opticlusterInformation;
    std::string headers;
};



#endif //OPTICLUSTERDATA_H
