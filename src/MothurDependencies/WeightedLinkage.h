//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef WEIGHTEDLINKAGE_H
#define WEIGHTEDLINKAGE_H
#include "Cluster.h"
#include <string>

class WeightedLinkage : public Cluster {
public:
    WeightedLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) override;
    std::string getTag() override;

private:
    unsigned long saveRow;
    unsigned long saveCol;
};

#endif //WEIGHTEDLINKAGE_H
