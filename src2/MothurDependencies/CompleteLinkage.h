//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef COMPLETELINKAGE_H
#define COMPLETELINKAGE_H
#include "Cluster.h"
#include <string>

class CompleteLinkage final : public Cluster {
public:
    CompleteLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, const std::string &, float);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) override;
    std::string getTag() override;

private:

};

#endif //COMPLETELINKAGE_H
