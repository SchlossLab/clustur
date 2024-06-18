//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef COMPLETELINKAGE_H
#define COMPLETELINKAGE_H
#include "Cluster.h"


class CompleteLinkage : public Cluster {
public:
    CompleteLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell);
    std::string getTag();

private:

};

#endif //COMPLETELINKAGE_H
