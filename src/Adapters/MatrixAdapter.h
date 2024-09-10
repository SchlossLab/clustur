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
        const std::vector<double> &dataValues, double cutOff, bool isSimularity, const CountTableAdapter& table);
    ~MatrixAdapter() = default;
    ReadPhylipMatrix* ReadPhylipFile(const std::string& path) const;
    SparseDistanceMatrix* CreateSparseMatrix();
    ListVector* GetListVector() const {return phylipReader->getListVector();}
    SparseDistanceMatrix* GetSpareDistanceMatrix() const {return phylipReader->getDMatrix();}

    bool CreatePhylipFile(const std::string &saveFileLocation);
    std::vector<RowData> DistanceMatrixToSquareMatrix();

private:

    ReadPhylipMatrix* phylipReader;
    SparseDistanceMatrix* spareDistanceMatrix = nullptr;
    CountTableAdapter countTable;
    std::vector<int> xPosition;
    std::vector<int> yPosition;
    std::vector<double> data;
    std::vector<std::string> matrixNames;
};



#endif //MATRIXADAPTER_H
