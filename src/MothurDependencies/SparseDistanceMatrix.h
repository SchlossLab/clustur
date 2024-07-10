//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef SPARSEDISTANCEMATRIX_H
#define SPARSEDISTANCEMATRIX_H

#include "PDistCell.h"

#include <vector>

#include "Utils.h"


class ListVector;



class SparseDistanceMatrix {

public:
    SparseDistanceMatrix();
    ~SparseDistanceMatrix(){ clear(); }
    int getNNodes();
    unsigned long getSmallestCell(unsigned long& index);		//Return the cell with the smallest distance
    float getSmallDist();
    static bool heapComparator(const PDistCell& a, const PDistCell& b);
    int rmCell(unsigned long, unsigned long);
    int updateCellCompliment(unsigned long, unsigned long);
    void resize(const unsigned long n) { seqVec.resize(n);  }
    void clear();
    void addCell(unsigned long, PDistCell);
    int addCellSorted(unsigned long, PDistCell);
    std::vector<std::vector<PDistCell> > seqVec;
    bool print();

private:
    PDistCell smallCell;				//The cell with the smallest distance
    int numNodes;

    bool sorted;
    int sortSeqVec();
    int sortSeqVec(int);
    float smallDist, aboveCutoff;

    Utils util;

};

#endif //SPARSEDISTANCEMATRIX_H
