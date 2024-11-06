//
// Created by Gregory Johnson on 6/28/24.
//

#ifndef MATRIXADAPTER_H
#define MATRIXADAPTER_H
#include <unordered_map>

#include "CountTableAdapter.h"
#include "../MothurDependencies/ReadPhylipMatrix.h"


class MatrixAdapter {
public:
    MatrixAdapter(const std::vector<int> &iIndexes, const std::vector<int> &jIndexes,
        const std::vector<double> &dataValues, double cutOff, bool isSimularity, CountTableAdapter table);
    ~MatrixAdapter() = default;
    bool CreatePhylipFile(const std::string &saveFileLocation);
    bool CreateColumnDataFile(const std::string &saveFileLocation);
    SparseDistanceMatrix DistanceMatrixToSparseMatrix();
    ListVector CreateListVector() const;

private:
    double cutoff;
    CountTableAdapter countTable;
    std::vector<int> xPosition;
    std::vector<int> yPosition;
    std::vector<double> data;
    std::vector<std::string> matrixNames;
};



#endif //MATRIXADAPTER_H
