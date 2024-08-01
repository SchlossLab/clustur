//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef SHAREDFILEBUILDER_H
#define SHAREDFILEBUILDER_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <set>
#include <list>

#include "ListVector.h"
#include "RSparseMatrix.h"
#include "SharedFile.h"
#include "OtuAbundancePair.h"

class SharedFileBuilder {
public:
    explicit SharedFileBuilder(const RSpraseMatrix& sparseMatrix);
    void CreateNameMap();
    void CreateCountTable();
    void OutputFiles() const;
    SharedFile* BuildSharedFile(ListVector*);

private:
    RSpraseMatrix matrix;
    std::unordered_map<int, std::set<int>> nameMap;
    std::unordered_map<int, int> countTable;
    Utils utils;
    float cutoff = 0.1;
};





#endif //SHAREDFILEBUILDER_H
