//
// Created by Gregory Johnson on 7/12/24.
//

#ifndef MATRIXADAPTERTESTFIXTURE_H
#define MATRIXADAPTERTESTFIXTURE_H
#include "TestFixture.h"
#include <vector>
#include "../Adapters/MatrixAdapter.h"


class MatrixAdapterTestFixture final: public TestFixture {
public:
    bool TestReadPhylipFile(bool expectedResult);
    bool TestCreateSparseMatrix(bool expectedResult);
    bool TestCreatePhylipFile(bool expectedResult);
    bool TestDistanceMatrixToSquareMatrix(int expectedResult);
    bool TestGetListVector(bool createSparseFirst, bool expectedResult);
    bool TestGetSparseDistanceMatrix(bool createSparseFirst, bool expectedResult);
    ~MatrixAdapterTestFixture() override;
private:

    void Setup() override;
    void TearDown() override;
    MatrixAdapter* matrixAdapter = nullptr;
};



#endif //MATRIXADAPTERTESTFIXTURE_H
