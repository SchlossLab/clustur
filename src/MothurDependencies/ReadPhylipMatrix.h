//
// Created by Gregory Johnson on 6/20/24.
//

#ifndef READPHYLIPMATRIX_H
#define READPHYLIPMATRIX_H

#include <fstream>
#include <string>
#include <vector>
#include "../RowData.h"
#include <unordered_set>
#include "ListVector.h"
#include "SparseDistanceMatrix.h"
#include "../Adapters/DistanceFileReader.h"


class ReadPhylipMatrix final : public DistanceFileReader {

public:
    ReadPhylipMatrix(double, bool);
    ReadPhylipMatrix() = default;
    // ~ReadPhylipMatrix() override = default;
    bool Read(const std::string&) override;
private:
    std::ifstream fileHandle;
};

#endif //READPHYLIPMATRIX_H
