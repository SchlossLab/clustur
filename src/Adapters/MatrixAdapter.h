//
// Created by Gregory Johnson on 6/28/24.
//

#ifndef MATRIXADAPTER_H
#define MATRIXADAPTER_H
#include "../MothurDependencies/ReadPhylipMatrix.h"


class MatrixAdapter {
public:
    MatrixAdapter();
    ~MatrixAdapter();
    void ReadPhylipFile();
    void ReadRowDataStructure();
    void CreatePhylipFile();

private:
    ReadPhylipMatrix* phylipMatrix;
};



#endif //MATRIXADAPTER_H
