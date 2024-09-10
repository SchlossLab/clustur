//
// Created by Gregory Johnson on 6/28/24.
//
#include "Adapters/MatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include <map>

MatrixAdapter::MatrixAdapter(const std::vector<int> &iIndexes, const std::vector<int> &jIndexes,
                             const std::vector<double> &dataValues, const double cutOff, const bool isSimularity,
                             const CountTableAdapter& table): countTable(table),xPosition(iIndexes), yPosition(jIndexes),
                            data(dataValues) {
    phylipReader = new ReadPhylipMatrix(cutOff, isSimularity);
}

ReadPhylipMatrix* MatrixAdapter::ReadPhylipFile(const std::string &path) const {
    if(path.empty())
        return nullptr;
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

bool MatrixAdapter::CreatePhylipFile(const std::string &saveFileLocation) {\
    if(saveFileLocation.empty())
        return false;
    const auto matrix = DistanceMatrixToSquareMatrix();
    int count = 0;
    const size_t size = matrixNames.size();
    std::string distanceString = "\t" + std::to_string(size) + "\n";
    for (const auto &cells: matrix) {
        distanceString += matrixNames[count];
        for (int i = 0; i < count; i++) {
            distanceString += "\t" + std::to_string(cells.rowValues[i]);
        }
        distanceString += "\n";
        count++;
    }
    std::ofstream writeOut(saveFileLocation);
    if (!writeOut.is_open()) {
        TestHelper::Print("Failed to open: \n");
    }
    writeOut << distanceString;
    writeOut.close();
    return true;
}

std::vector<RowData> MatrixAdapter::DistanceMatrixToSquareMatrix() {
    // The indexes are +1, i need to push them back so that 1 -> 0, 2-> 1, etc (name map maybe?)
    std::set<std::string> names;
    const size_t nSeqs = data.size();
    std::map<int, RowData> dataList;
    std::unordered_map<int, int> positionsOfIndexs;
    std::unordered_map<int, std::string> positionsToNames;
    for (int i = 0; i < nSeqs; i++) {
        names.insert(countTable.GetNameByIndex(i));
        positionsToNames[xPosition[i]] = countTable.GetNameByIndex(i); // Not going to work, I need a way to link my names to the sparse matix indices
    }

    const size_t nameSize = names.size();
    matrixNames = std::vector<std::string>(nameSize);
    for (int i = 0; i < nameSize; i++) {
        positionsOfIndexs[xPosition[i]] = i;
        matrixNames[i] = positionsToNames[xPosition[i]];
        dataList[i].name = positionsToNames[xPosition[i]];
        dataList[i].rowValues = std::vector<double>(nameSize, 0);
    }

    for (size_t i = 0; i < nSeqs;  i++) {
        int xIndex = positionsOfIndexs[xPosition[i]]; // Coming from r -> c++, indeces start at 1 in r
        int yIndex = positionsOfIndexs[yPosition[i]];
        if(data[i] < 0) {
            data[i] = 0;
        }
        const double currentValueX = dataList[yIndex].rowValues[xIndex];
        if(currentValueX != 0){ // We already set the value and this is a sparse matrix. 
            continue;           // WE do not need to reset the values back to zero.
        }                       // This is a catch all in the case of a sparse and square matrix
        // Since the indexes were reverting back to zero, if the values were found again,
        // like 2,4 = 0.3, but 4,2 = 0 was found, (its a sparse matrix) so we do not change back the value.
        dataList[xIndex].rowValues[yIndex] = data[i];
        dataList[yIndex].rowValues[xIndex] = data[i];
    }
    std::vector<RowData> sequenceData(names.size());
    int index = 0;
    for (const auto &row: dataList) {
        sequenceData[index++] = row.second;
    }
    return sequenceData;
}
