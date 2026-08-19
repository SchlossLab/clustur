//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_SPARSEDISTANCEMATRIX_H
#define REFACTOR_SPARSEDISTANCEMATRIX_H
#define MOTHURMAX 1e6


#include <vector>
#include "../MothurDependencies/Utils.h"
#include "PDistCell.h"

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
    SparseDistanceMatrix(size_t size);

    ~SparseDistanceMatrix(){ clear(); }
    [[nodiscard]] int getNNodes() const;
    unsigned long getSmallestCell(unsigned long& row);		//Return the cell with the smallest distance
    [[nodiscard]] float getSmallDist() const;
    static bool heapComparator(const PDistCell& a, const PDistCell& b);
    int rmCell(unsigned long, unsigned long);
    void updateCellCompliment(unsigned long, unsigned long);
    void resize(const unsigned long n) {
        if(n > seqVec.max_size() || n <= 0)
            return;
        seqVec.resize(n);
    }
    void clear();
    void addCells(const SparseDistanceMatrix &other);
    void addCell(unsigned long, PDistCell);
    int addCellSorted(unsigned long, PDistCell);
    std::vector<std::vector<PDistCell> > seqVec;
    [[nodiscard]] bool print() const;
    void FilterSparseMatrix(float cutoff);

private:
    PDistCell smallCell;				//The cell with the smallest distance
    int numNodes;

    bool sorted;
    int sortSeqVec();
    int sortSeqVec(unsigned long);
    float smallDist, aboveCutoff;

    Utils util;

};

#endif //REFACTOR_SPARSEDISTANCEMATRIX_H