//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef SINGLELINKAGE_H
#define SINGLELINKAGE_H
#include "Cluster.h"
#include <string>

#include "ClusterMethod.h"

class SingleLinkage final : public Cluster, public ClusterMethod {
public:
    SingleLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, const std::string&, float);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) override;
    std::string getTag() override;

    ClusterExport * Execute() override;
};




#endif //SINGLELINKAGE_H
