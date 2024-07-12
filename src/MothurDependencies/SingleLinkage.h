//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef SINGLELINKAGE_H
#define SINGLELINKAGE_H
#include "Cluster.h"
#include <string>

class SingleLinkage : public Cluster {
public:
    SingleLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    //void update(double&);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) override;
    std::string getTag() override;

private:

};




#endif //SINGLELINKAGE_H
