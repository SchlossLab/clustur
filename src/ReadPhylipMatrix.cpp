//
// Created by Gregory Johnson on 6/20/24.
//

#include "MothurDependencies/ReadPhylipMatrix.h"
#include "MothurDependencies/ListVector.h"
#include "TestHelpers/TestHelper.h"

ReadPhylipMatrix::ReadPhylipMatrix(const double cutoff, const bool simularityMatrix):cutoff(cutoff), sim(simularityMatrix) {
    DMatrix = new SparseDistanceMatrix();
    list = new ListVector();
}
bool ReadPhylipMatrix::read(const std::string& filePath) {

    fileHandle.open(filePath);
    if(!fileHandle.is_open())
        return false;

    float distance;
    int square = 0;
    std::string name;
    std::vector<std::string> matrixNames;


    std::string numTest;
    fileHandle >> numTest >> name;
    const int nseqs = std::stoi(numTest);
    matrixNames.push_back(name);


    list = new ListVector(nseqs);
    list->set(0, name);


    char d;
    while ((d = fileHandle.get()) != EOF) {
        if (isalnum(d)) {
            square = 1;
            fileHandle.putback(d);
            for (int i = 0; i < nseqs; i++) {
                fileHandle >> distance;
            }
            break;
        }
        if (d == '\n') {
            square = 0;
            break;
        }
    }

    DMatrix->resize(nseqs);
    if (square == 0) {

        for (int i = 1; i < nseqs; i++) {
            fileHandle >> name;
            matrixNames.push_back(name);

            list->set(i, name);
            for (int j = 0; j < i; j++) {
                fileHandle >> distance;

                if (util.isEqual(distance, -1)) { distance = 1000000; } else if (sim) { distance = 1 - distance; }
                //user has entered a sim matrix that we need to convert.

                if (distance <= cutoff) {
                    const PDistCell value(i, distance);
                    DMatrix->addCell(j, value);
                }
            }

        }
    }
    else {

        for (int i = 1; i < nseqs; i++) {
            fileHandle >> name;
            matrixNames.push_back(name);

            list->set(i, name);
            for (int j = 0; j < nseqs; j++) {
                fileHandle >> distance;


                if (util.isEqual(distance, -1)) { distance = 1000000; } else if (sim) { distance = 1 - distance; }
                //user has entered a sim matrix that we need to convert.

                if (distance <= cutoff && j < i) {
                    const PDistCell value(i, distance);
                    DMatrix->addCell(j, value);
                }
            }

        }
    }
    list->setLabel("0");
    fileHandle.close();

    return true;
}

std::vector<RowData> ReadPhylipMatrix::readToRowData(const std::string& filePath) {

    fileHandle.open(filePath);
    if(!fileHandle.is_open())
        return {};

    float distance;
    std::string name;
    std::vector<std::string> matrixNames;


    std::string numTest;
    fileHandle >> numTest >> name;
    const int nseqs = std::stoi(numTest);
    matrixNames.push_back(name);
    std::vector<RowData> phylipRowData(nseqs);
    phylipRowData[0].name = name;
    list = new ListVector(nseqs);
    list->set(0, name);

    for (int i = 0; i < nseqs; i++) {
        phylipRowData[i].rowValues = std::vector<double>(nseqs);
    }

    for (int i = 1; i < nseqs; i++) {
        fileHandle >> name;
        matrixNames.push_back(name);
        phylipRowData[i].name = name;
        list->set(i, name);
        for (int j = 0; j < i; j++) {
            fileHandle >> distance;
            phylipRowData[i].rowValues[j] = distance;
            phylipRowData[j].rowValues[i] = distance;
        }
    }

    fileHandle.close();

    return phylipRowData;
}


bool ReadPhylipMatrix::read(const std::vector<RowData> &rowData) {
    if (rowData.empty())
        return false;
    std::string name = rowData[0].name;
    std::vector<std::string> matrixNames;
    const int nseqs = static_cast<int>(rowData.size());
    list = new ListVector(nseqs);
    list->set(0, name);

    DMatrix->resize(nseqs);
    for (int i = 1; i < nseqs; i++) {
        name = rowData[i].name;
        list->set(i, name);
        matrixNames.push_back(name); // Square matrix respresented as a sparse matrix. This causes a myraid of problems.
        // Where it should be giving me the distance between 1 and 2, its not giving me any distance, I removed them.
        for (int j = 0; j < i; j++) {
            auto distance = static_cast<float>(rowData[i].rowValues[j]);
            const bool equalivance = util.isEqual(distance, -1);
            if (equalivance) {
                distance = 1000000;
            } else if (sim) {
                distance = 1.0f - distance;
            }

            if (distance <= cutoff) {
                const PDistCell value(i, distance);
                DMatrix->addCell(j, value);
            }
        }
    }
    list->setLabel("0");
    return true;
}
