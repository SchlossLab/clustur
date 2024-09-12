//
// Created by Gregory Johnson on 6/20/24.
//

#ifndef READPHYLIPMATRIX_H
#define READPHYLIPMATRIX_H

#include <fstream>
#include <string>
#include <vector>
#include "../RowData.h"
#include "ListVector.h"
#include "SparseDistanceMatrix.h"


class ReadPhylipMatrix  {

public:
    ReadPhylipMatrix(double, bool);
    ReadPhylipMatrix() = default;
    // ReadPhylipMatrix(std::string, bool);
    ~ReadPhylipMatrix() {};
    SparseDistanceMatrix* getDMatrix() const { return DMatrix;}
    ListVector* getListVector()	const {	return list;}
    bool read(const std::string&);
    bool read(const std::vector<RowData>& rowData);
    std::vector<RowData> readToRowData(const std::string&);
private:
    std::ifstream fileHandle;
    SparseDistanceMatrix* DMatrix{};
    ListVector* list{};
    double cutoff = 0;
    bool sim = true;
    Utils util;
};

#endif //READPHYLIPMATRIX_H
