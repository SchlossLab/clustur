//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef AVERAGELINKAGE_H
#define AVERAGELINKAGE_H
#include "Cluster.h"
#include <string>

class AverageLinkage final: public Cluster {
public:
    AverageLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) override;
   std::string getTag() override;

private:
    unsigned long saveRow;
    unsigned long saveCol;
    unsigned long rowBin;
    unsigned long colBin;
    unsigned long totalBin;

};


#endif //AVERAGELINKAGE_H
