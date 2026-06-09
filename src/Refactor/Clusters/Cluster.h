//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_CLUSTER_H
#define REFACTOR_CLUSTER_H




#include "../DataStructures/SparseDistanceMatrix.h"
#include "../DataStructures/PDistCell.h"
#include "../DataStructures/RAbundVector.h"
#include "../DataStructures/ListVector.h"

// typedef unsigned long ull;
class Cluster {

public:
    Cluster(const RAbundVector&, const ListVector&, const SparseDistanceMatrix&, float, std::string, float);
    Cluster() = default;
    virtual ~Cluster() = default;
    virtual bool update(double&);
    virtual std::string getTag() = 0;
    virtual std::map<std::string, int> getSeqtoBin()  {  return seq2Bin;	}
    virtual bool updateDistance(PDistCell& colCell, PDistCell& rowCell) = 0;
protected:
    virtual bool clusterBins();
    virtual bool clusterNames();
    virtual bool updateMap();


    RAbundVector rabund{};
    ListVector list{};
    SparseDistanceMatrix dMatrix{};

    unsigned long smallRow{}, smallCol{}, nRowCells{}, nColCells{};
    float cutoff{};
    std::string method;
    float smallDist{}, adjust{};
    bool mapWanted{};
    std::map<std::string, int> seq2Bin;

    Utils util;
};


#endif //REFACTOR_CLUSTER_H