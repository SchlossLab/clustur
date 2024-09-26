//
// Created by Gregory Johnson on 9/13/24.
//

#ifndef COLUMNDISTANCEMATRIXREADER_H
#define COLUMNDISTANCEMATRIXREADER_H
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include "ListVector.h"
#include "SparseDistanceMatrix.h"
#include "../RowData.h"
#include "../Adapters/CountTableAdapter.h"

class ColumnDistanceMatrixReader final {
public:
    ColumnDistanceMatrixReader(double cutoff, bool isSimularity);
    bool Read(const CountTableAdapter &countTable, const std::string &filePath);
    std::vector<RowData> readToRowData(const CountTableAdapter &countTable, const std::string& filePath);
    ListVector* GetListVector() const {return list;}
    SparseDistanceMatrix* GetSparseDataMatrix() const {return sparseMatrix;}

private:
    SparseDistanceMatrix* sparseMatrix{};
    ListVector* list{};
    bool sim;
    double cutoff;
};



#endif //COLUMNDISTANCEMATRIXREADER_H
