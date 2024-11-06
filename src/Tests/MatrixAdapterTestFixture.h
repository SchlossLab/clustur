//
// Created by Gregory Johnson on 7/12/24.
//

#ifndef MATRIXADAPTERTESTFIXTURE_H
#define MATRIXADAPTERTESTFIXTURE_H
#include "TestFixture.h"
#include <vector>
#include <Rcpp.h>
#include "../Adapters/MatrixAdapter.h"


class MatrixAdapterTestFixture final: public TestFixture {
public:
    bool TestCreateSparseMatrix(size_t expectedResult);
    bool TestCreatePhylipFile(bool expectedResult);
    bool TestGetListVector(bool createSparseFirst, bool expectedResult);
    ~MatrixAdapterTestFixture() override;
private:

    void Setup() override;
    void TearDown() override;
    MatrixAdapter* matrixAdapter = nullptr;
};



#endif //MATRIXADAPTERTESTFIXTURE_H
