//
// Created by Gregory Johnson on 6/20/24.
//

#include "MothurDependencies/ReadPhylipMatrix.h"
#include "MothurDependencies/ListVector.h"

ReadPhylipMatrix::ReadPhylipMatrix(const double cutoff, const bool simularityMatrix):DistanceFileReader() {
    this->cutoff = cutoff;
    this->sim = simularityMatrix;
    sparseMatrix = new SparseDistanceMatrix();
    list = new ListVector();
}

bool ReadPhylipMatrix::Read(const std::string& filePath) {

    fileHandle.open(filePath);
    if(!fileHandle.is_open())
        return false;
    Utils util;
    float distance;
    int square = 0;
    std::string name;
    const Utils utils;
    const std::vector<std::string> samples = countTable.GetSamples();
    const std::unordered_set<std::string> sampleContainer(samples.begin(), samples.end());


    std::string numTest;
    fileHandle >> numTest >> name;
    const int nseqs = std::stoi(numTest);


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

    sparseMatrix->resize(nseqs);
    if (square == 0) {

        for (int i = 1; i < nseqs; i++) {
            fileHandle >> name;
            if(sampleContainer.find(name) == sampleContainer.end())
                utils.CheckForDistanceFileError({name});

            list->set(i, name);
            for (int j = 0; j < i; j++) {
                fileHandle >> distance;

                if (util.isEqual(distance, -1)) { distance = 1000000; } else if (sim) { distance = 1 - distance; }
                //user has entered a sim matrix that we need to convert.

                if (distance <= cutoff) {
                    const PDistCell value(i, distance);
                    sparseMatrix->addCell(j, value);
                }
            }

        }
    }
    else {

        for (int i = 1; i < nseqs; i++) {
            fileHandle >> name;
            if(sampleContainer.find(name) == sampleContainer.end())
                utils.CheckForDistanceFileError({name});

            list->set(i, name);
            for (int j = 0; j < nseqs; j++) {
                fileHandle >> distance;


                if (util.isEqual(distance, -1)) { distance = 1000000; } else if (sim) { distance = 1 - distance; }
                //user has entered a sim matrix that we need to convert.

                if (distance <= cutoff && j < i) {
                    const PDistCell value(i, distance);
                    sparseMatrix->addCell(j, value);
                }
            }

        }
    }
    list->setLabel("0");
    fileHandle.close();

    return true;
}
