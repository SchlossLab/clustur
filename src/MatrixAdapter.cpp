//
// Created by Gregory Johnson on 6/28/24.
//
#include "Adapters/MatrixAdapter.h"


MatrixAdapter::MatrixAdapter(const std::vector<int> &iIndexes, const std::vector<int> &jIndexes,
                             const std::vector<double> &dataValues, const double cutoff) {
    xPosition = iIndexes;
    yPosition = jIndexes;
    data = dataValues;
    this->cutoff = cutoff;
    phylipReader = new ReadPhylipMatrix(cutoff);
}

ReadPhylipMatrix* MatrixAdapter::ReadPhylipFile(const std::string &path) {
    phylipReader->read(path);
    return phylipReader;
}

SparseDistanceMatrix* MatrixAdapter::CreateSparseMatrix() {
    if(spareDistanceMatrix != nullptr)
        return spareDistanceMatrix;
    const auto phylipMatrix = DistanceMatrixToSquareMatrix();
    phylipReader->read(phylipMatrix);
    spareDistanceMatrix = phylipReader->getDMatrix();

    return spareDistanceMatrix;
}

void MatrixAdapter::CreatePhylipFile(const std::string& saveFileLocation) {
    auto matrix = CreateSparseMatrix();
    int count = 0;
    const size_t size = matrixNames.size();
    std::string distanceString = "\t" + std::to_string(size) + "\n";
    for (const auto &cells: matrix->seqVec) {
        if (cells.empty()) continue;
        distanceString += matrixNames[count++];
        int indexCounter = 0;
        for (int i = 0; i < size; i++) {
            if (!cells.empty() && cells[indexCounter].index == i) {
                distanceString += "\t" + std::to_string(cells[indexCounter++].dist);
                continue;
            }
            distanceString += "\t0";
        }
        distanceString += "\n";
    }
    std::ofstream writeOut(saveFileLocation);
    if (!writeOut.is_open()) {
        std::cout << "Failed to open: \n" << std::endl;
    }
    writeOut << distanceString;
    writeOut.close();
    delete(matrix);

}

std::vector<RowData> MatrixAdapter::DistanceMatrixToSquareMatrix() {
    // The indexes are +1, i need to push them back so that 1 -> 0, 2-> 1, etc (name map maybe?)
    std::set<int> names;
    const size_t nSeqs = data.size();
    std::unordered_map<int, RowData> dataList;
    std::unordered_map<int, int> positionsOfIndexs;
    for (size_t i = 0; i < nSeqs; i++) {
        names.insert(xPosition[i]);
        names.insert(yPosition[i]);
    }

    auto nameIter = names.begin();
    const size_t nameSize = names.size();
    matrixNames = std::vector<std::string>(nameSize);
    for (int i = 0; i < nameSize; i++) {
        positionsOfIndexs[*nameIter] = i;
        matrixNames[i] = std::to_string(*nameIter);
        dataList[i].name = std::to_string(*nameIter++);
        dataList[i].rowValues = std::vector<double>(nameSize);
    }

    for (size_t i = 0; i < nSeqs;  i++) {
        int xIndex = positionsOfIndexs[xPosition[i]];
        int yIndex = positionsOfIndexs[yPosition[i]];


        if(data[i] < 0) {
            data[i] = 0;
        }
        dataList[xIndex].rowValues[yIndex] = data[i];
        dataList[yIndex].rowValues[xIndex] = data[i];
    }
    std::vector<RowData> sequenceData(names.size());
    int index = 0;
    for (const auto &row: dataList) {
        sequenceData[index++] = row.second;
    }
    std::sort(sequenceData.begin(), sequenceData.end(), RowData::SortComparison);
    return sequenceData;
}
