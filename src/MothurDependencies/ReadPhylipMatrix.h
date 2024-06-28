//
// Created by Gregory Johnson on 6/20/24.
//

#ifndef READPHYLIPMATRIX_H
#define READPHYLIPMATRIX_H
#include <fstream>

#include "../RowData.h"
#include "ListVector.h"
#include "sparsedistancematrix.h"


class ReadPhylipMatrix  {

public:
    ReadPhylipMatrix(double);
    // ReadPhylipMatrix(std::string, bool);
    ~ReadPhylipMatrix() {};
    SparseDistanceMatrix* getDMatrix() const { return DMatrix;}
    ListVector* getListVector()	const {	return list;}
    int read(const std::string&);
    int read(const std::vector<RowData>& rowData);
private:
    std::ifstream fileHandle;
    SparseDistanceMatrix* DMatrix;
    ListVector* list;
    double cutoff;
    bool sim = false;
    Utils util;
};

#endif //READPHYLIPMATRIX_H
