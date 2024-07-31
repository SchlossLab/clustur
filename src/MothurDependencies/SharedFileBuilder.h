//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef SHAREDFILEBUILDER_H
#define SHAREDFILEBUILDER_H

#include <vector>
#include <unordered_map>

#include "RSparseMatrix.h"


class SharedFileBuilder {
public:
    explicit SharedFileBuilder(const RSpraseMatrix& sparseMatrix);
    void CreateNameMap();
    void CreateCountTable();
    void BuildSharedFile();
private:
    RSpraseMatrix matrix;
    float cutoff = 0.1;
};



#endif //SHAREDFILEBUILDER_H
