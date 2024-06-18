//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef SINGLELINKAGE_H
#define SINGLELINKAGE_H
#include "Cluster.h"


class SingleLinkage : public Cluster {
public:
    SingleLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    //void update(double&);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell);
    std::string getTag();

private:

};




#endif //SINGLELINKAGE_H
