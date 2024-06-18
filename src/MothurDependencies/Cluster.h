//
// Created by Gregory Johnson on 6/14/24.
//

#ifndef CLUSTER_H
#define CLUSTER_H



#include "sparsedistancematrix.h"
#include "optimatrix.h"
#include "PDistCell.h"
#include "RAbundVector.h"

class ListVector;

// typedef unsigned long ull;
class Cluster {

public:
	Cluster(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    Cluster() {}
    virtual ~Cluster() = default;
    virtual bool update(double&);
	virtual std::string getTag() = 0;
	virtual void setMapWanted(bool m);
	virtual std::map<std::string, int> getSeqtoBin()  {  return seq2Bin;	}
protected:


	virtual bool updateDistance(PDistCell& colCell, PDistCell& rowCell) = 0;

	virtual void clusterBins();
	virtual void clusterNames();
	virtual void updateMap();


	RAbundVector* rabund;
	ListVector* list;
	SparseDistanceMatrix* dMatrix;

	unsigned long smallRow, smallCol, nRowCells, nColCells;
	float smallDist, adjust, cutoff;
	bool mapWanted;
	std::map<std::string, int> seq2Bin;
	std::string method;

    Utils util;
};


#endif //CLUSTER_H
