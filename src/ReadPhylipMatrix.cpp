//
// Created by Gregory Johnson on 6/20/24.
//

#include "MothurDependencies/ReadPhylipMatrix.h"
#include "MothurDependencies/ListVector.h"

ReadPhylipMatrix::ReadPhylipMatrix(const double cutoff, const bool simularityMatrix) {
    DMatrix = new SparseDistanceMatrix();
    list = new ListVector();
    this->cutoff = cutoff;
    sim = simularityMatrix;
}

bool ReadPhylipMatrix::read(const std::string& filePath) {

    fileHandle.open(filePath);
    if(!fileHandle.is_open())
        return false;

    const ListVector *nameMap = nullptr;
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
        int index = 0;

        for (int i = 1; i < nseqs; i++) {
            fileHandle >> name;
            matrixNames.push_back(name);

            list->set(i, name);
            for (int j = 0; j < i; j++) {
                fileHandle >> distance;

                if (util.isEqual(distance, -1)) { distance = 1000000; } else if (sim) { distance = 1.0 - distance; }
                //user has entered a sim matrix that we need to convert.

                if (distance <= cutoff) {
                    const PDistCell value(i, distance);
                    DMatrix->addCell(j, value);
                }
                index++;
            }

        }
    }
    else {
        int index = nseqs;

        for (int i = 1; i < nseqs; i++) {
            fileHandle >> name;
            matrixNames.push_back(name);

            list->set(i, name);
            for (int j = 0; j < nseqs; j++) {
                fileHandle >> distance;


                if (util.isEqual(distance, -1)) { distance = 1000000; } else if (sim) { distance = 1.0 - distance; }
                //user has entered a sim matrix that we need to convert.

                if (distance <= cutoff && j < i) {
                    const PDistCell value(i, distance);
                    DMatrix->addCell(j, value);
                }
                index++;
            }

        }
    }
    list->setLabel("0");
    fileHandle.close();

    return true;
}

bool ReadPhylipMatrix::read(const std::vector<RowData> &rowData) {
    if (rowData.empty())
        return false;
    std::string name = rowData[0].name;
    std::vector<std::string> matrixNames;
    const size_t nseqs = rowData.size();
    list = new ListVector(nseqs);
    list->set(0, name);

    DMatrix->resize(nseqs);
    for (int i = 1; i < nseqs; i++) {
        name = rowData[i].name;
        list->set(i, name);
        matrixNames.push_back(name);
        for (int j = 0; j < nseqs; j++) {
            float distance = rowData[i].rowValues[j];
            const bool equalivance = util.isEqual(distance, -1);
            if (equalivance) {
                distance = 1000000;
            } else if (sim) {
                distance = 1.0f - distance;
            }
            if (distance <= cutoff && j < i) {
                const PDistCell value(i, distance);
                DMatrix->addCell(j, value);
            }
        }
    }
    list->setLabel("0");
    return true;
}
