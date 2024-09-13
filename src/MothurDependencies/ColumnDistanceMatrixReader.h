//
// Created by Gregory Johnson on 9/13/24.
//

#ifndef COLUMNDISTANCEMATRIXREADER_H
#define COLUMNDISTANCEMATRIXREADER_H
#include <string>
#include <fstream>
#include <map>
#include "ListVector.h"
#include "SparseDistanceMatrix.h"
#include "../RowData.h"
#include "../Adapters/CountTableAdapter.h"

class ColumnDistanceMatrixReader final {
public:
    ColumnDistanceMatrixReader(std::string  path, double cutoff, bool isSimularity);
    bool Read(const CountTableAdapter& countTable);
    std::vector<RowData> readToRowData(const CountTableAdapter &countTable);
    ListVector* GetListVector() const {return list;}
    SparseDistanceMatrix* GetSparseDataMatrix() const {return sparseMatrix;}

private:
    SparseDistanceMatrix* sparseMatrix{};
    ListVector* list{};
    std::string filePath;
    bool sim;
    double cutoff;
};



#endif //COLUMNDISTANCEMATRIXREADER_H
