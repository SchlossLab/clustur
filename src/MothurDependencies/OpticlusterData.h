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

class OpticlusterData : public ClusterExport{
public:
    explicit OpticlusterData(std::string header, const double cutoff):headers(std::move(header)){}
    void AddToData(OptiClusterInformation&);
    std::string Print() const override;
    std::string GetLabel(int) const override;

private:
    std::vector<OptiClusterInformation> opticlusterInformation;
    std::string headers;
};



#endif //OPTICLUSTERDATA_H
