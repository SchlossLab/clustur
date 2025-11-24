//
// Created by Gregory Johnson on 7/12/24.
//

#include "Tests/MatrixAdapterTestFixture.h"

bool MatrixAdapterTestFixture::TestCreateSparseMatrix(const MatrixAdapter& adapter, const size_t expectedResult) {
    Setup();
    const auto dat = adapter.CreateSparseMatrix();
    const auto result = dat.seqVec.size();
    TearDown();
    return result == expectedResult;
}

bool MatrixAdapterTestFixture::TestGetListVector(const MatrixAdapter& adapter, const bool expectedResult) {
    Setup();
    const auto result = adapter.CreateListVector().size() > 0;
    TearDown();
    return result == expectedResult;
}

MatrixAdapterTestFixture::~MatrixAdapterTestFixture() {
    TearDown();
}

void MatrixAdapterTestFixture::Setup() {
    // Do nothing
}

void MatrixAdapterTestFixture::TearDown() {
    // Do nothing
}
