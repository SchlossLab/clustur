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
#include "../Adapters/DistanceFileReader.h"


class ReadPhylipMatrix final : public DistanceFileReader {

public:
    ReadPhylipMatrix(double, bool);
    ReadPhylipMatrix() = default;
    ~ReadPhylipMatrix() override = default;
    SparseDistanceMatrix* GetSparseMatrix() const override { return DMatrix;}
    ListVector* GetListVector()	const override{	return list;}
    bool read(const std::string&);
    bool read(const std::vector<RowData>& rowData);
    std::vector<RowData> ReadToRowData(const std::string&) override;

    std::vector<RowData> ReadToRowData(const CountTableAdapter &adapter, const std::string &filePath) override { return {};};

private:
    std::ifstream fileHandle;
    SparseDistanceMatrix* DMatrix{};
    ListVector* list{};
    double cutoff = 0;
    bool sim = true;
    Utils util;
};

#endif //READPHYLIPMATRIX_H
