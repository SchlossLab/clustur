//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef AVERAGELINKAGE_H
#define AVERAGELINKAGE_H
#include "Cluster.h"


class AverageLinkage : public Cluster {
public:
    AverageLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell);
   std::string getTag();

private:
    int saveRow;
    int saveCol;
    int rowBin;
    int colBin;
    int totalBin;

};


#endif //AVERAGELINKAGE_H
