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
    ReadPhylipMatrix(std::string, float);
    //ReadPhylipMatrix(std::string, bool);
    ~ReadPhylipMatrix() {};
    SparseDistanceMatrix* getDMatrix() const { return DMatrix;}
    ListVector* getListVector()	const {	return list;}
    int read();
    int read(const std::vector<RowData>& rowData);
private:
    std::ifstream fileHandle;
    std::string distFile;
    SparseDistanceMatrix* DMatrix;
    ListVector* list;
    float cutoff;
    bool sim = false;
    Utils util;
};

#endif //READPHYLIPMATRIX_H
