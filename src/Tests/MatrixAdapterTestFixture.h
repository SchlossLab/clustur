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
    bool TestCreateSparseMatrix(const MatrixAdapter& adapter, size_t expectedResult);
    bool TestGetListVector(const MatrixAdapter& adapter, bool expectedResult);
    ~MatrixAdapterTestFixture() override;
private:

    void Setup() override;
    void TearDown() override;
};



#endif //MATRIXADAPTERTESTFIXTURE_H
