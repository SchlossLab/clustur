//
// Created by Gregory Johnson on 6/28/24.
//
#include "Adapters/MatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include <map>
#include <utility>

MatrixAdapter::MatrixAdapter(const std::vector<int> &iIndexes, const std::vector<int> &jIndexes,
                             const std::vector<double> &dataValues, const double cutOff, const bool isSimularity,
                             CountTableAdapter table): cutoff(cutOff), countTable(std::move(table)),xPosition(iIndexes), yPosition(jIndexes),
                            data(dataValues) {
}


bool MatrixAdapter::CreatePhylipFile(const std::string &saveFileLocation) {
    // if(saveFileLocation.empty())
    //     return false;
    // const auto matrix = DistanceMatrixToSquareMatrix();
    // int count = 0;
    // const size_t size = matrixNames.size();
    // std::string distanceString = "\t" + std::to_string(size) + "\n";
    // for (const auto &cells: matrix) {
    //     distanceString += matrixNames[count];
    //     for (int i = 0; i < count; i++) {
    //         distanceString += "\t" + std::to_string(cells.rowValues[i]);
    //     }
    //     distanceString += "\n";
    //     count++;
    // }
    // std::ofstream writeOut(saveFileLocation);
    // if (!writeOut.is_open()) {
    //     TestHelper::Print("Failed to open: \n");
    // }
    // writeOut << distanceString;
    // writeOut.close();
    return true;
}

bool MatrixAdapter::CreateColumnDataFile(const std::string &saveFileLocation) {
    if(saveFileLocation.empty())
        return false;
    const auto matrix = CreateSparseMatrix();
    std::string data;
    for (size_t i = 0; i < matrix.seqVec.size(); i++) {
        std::string firstCellName = matrixNames[i];
        for(size_t j = 0; j < matrixNames[i].size(); j++) {
            std::string otherCell = matrixNames[matrix.seqVec[i][j].index];
            data += firstCellName + "\t" + otherCell + "\t" + std::to_string(matrix.seqVec[i][j].dist) + "\n";
        }
    }
    std::ofstream writeOut(saveFileLocation);
    if (!writeOut.is_open()) {
        TestHelper::Print("Failed to open: \n");
    }
    writeOut << data;
    writeOut.close();
    return true;
}

SparseDistanceMatrix MatrixAdapter::CreateSparseMatrix() {
    // The indexes are +1, i need to push them back so that 1 -> 0, 2-> 1, etc (name map maybe?)
    std::set<std::string> names;
    SparseDistanceMatrix sparseMatrix;
    const int nSeqs = static_cast<int>(data.size());
    if(nSeqs <= 0)
        return {};
    std::unordered_map<int, int> positionsOfIndexs;
    std::unordered_map<int, std::string> positionsToNames;
    auto samples = countTable.GetSamples();
    names.insert(samples.begin(), samples.end());
    const int nameSize = static_cast<int>(names.size());
    sparseMatrix.resize(nameSize);
    if(static_cast<int>(xPosition.size()) > nameSize) { // There are values that should exist
        std::set<std::string> unknownNames;
        for(int i = nameSize; i < static_cast<int>(xPosition.size()); i++) {
            if(i >= nameSize + 2)
                break;
            unknownNames.insert(std::to_string(xPosition[i]));
        }
        const Utils util;
        util.CheckForDistanceFileError(unknownNames);
    }

    for (int i = 0; i < nameSize; i++) {
        positionsToNames[xPosition[i]] = countTable.GetNameByIndex(i); // Not going to work, I need a way to link my names to the sparse matix indices
    }


    matrixNames = std::vector<std::string>(nameSize);
    for (int i = 0; i < nameSize; i++) {
        positionsOfIndexs[xPosition[i]] = i;
        matrixNames[i] = positionsToNames[xPosition[i]];
    }


    for (int i = 0; i < nSeqs;  i++) {

        double currentDist = data[i];
        if(currentDist > cutoff) continue;
        if(currentDist < 0) {
            currentDist = 0;
        }
        const int xIndex = positionsOfIndexs[xPosition[i]]; // Coming from r -> c++, indeces start at 1 in r
        const int yIndex = positionsOfIndexs[yPosition[i]];

        // const double currentValueX = dataList[yIndex].rowValues[xIndex];
        // if(currentValueX != 0){ // We already set the value and this is a sparse matrix.
        //     continue;           // WE do not need to reset the values back to zero.
        // }                       // This is a catch all in the case of a sparse and square matrix
        // Since the indexes were reverting back to zero, if the values were found again,
        // like 2,4 = 0.3, but 4,2 = 0 was found, (its a sparse matrix) so we do not change back the value.
        if(xIndex > yIndex)
            sparseMatrix.addCell(yIndex, PDistCell(xIndex, static_cast<float>(currentDist)));
        else
            sparseMatrix.addCell(xIndex, PDistCell(yIndex, static_cast<float>(currentDist)));
    }
    return sparseMatrix;
}

ListVector MatrixAdapter::CreateListVector() const {
    ListVector vector;
    vector.resize(static_cast<int>(matrixNames.size()));
    int count = 0;
    for(const auto& name : matrixNames) {
        vector.set(count++, name);
    }
    return vector;
}
