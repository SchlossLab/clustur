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
    SparseDistanceMatrix(const SparseDistanceMatrix &other)
        : seqVec(other.seqVec),
          smallCell(other.smallCell),
          numNodes(other.numNodes),
          sorted(other.sorted),
          smallDist(other.smallDist),
          aboveCutoff(other.aboveCutoff),
          util(other.util) {}
    SparseDistanceMatrix();
    ~SparseDistanceMatrix(){ clear(); }
    int getNNodes() const;
    unsigned long getSmallestCell(unsigned long& row);		//Return the cell with the smallest distance
    float getSmallDist() const;
    static bool heapComparator(const PDistCell& a, const PDistCell& b);
    int rmCell(unsigned long, unsigned long);
    int updateCellCompliment(unsigned long, unsigned long);
    void resize(const unsigned long n) {
        if(n > seqVec.max_size() || n <= 0)
            return;
        seqVec.resize(n);
    }
    void clear();
    void addCell(unsigned long, PDistCell);
    int addCellSorted(unsigned long, PDistCell);
    std::vector<std::vector<PDistCell> > seqVec;
    bool print() const;

private:
    PDistCell smallCell;				//The cell with the smallest distance
    int numNodes;

    bool sorted;
    int sortSeqVec();
    int sortSeqVec(unsigned long);
    float smallDist, aboveCutoff;

    Utils util;

};

#endif //SPARSEDISTANCEMATRIX_H
