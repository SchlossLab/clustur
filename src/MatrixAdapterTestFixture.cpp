//
// Created by Gregory Johnson on 7/12/24.
//

#include "Tests/MatrixAdapterTestFixture.h"

bool MatrixAdapterTestFixture::TestReadPhylipFile(const bool expectedResult) {
    Setup();
    const bool result = matrixAdapter->ReadPhylipFile("") != nullptr;
    TearDown();
    return expectedResult == result;
}

bool MatrixAdapterTestFixture::TestCreateSparseMatrix(const bool expectedResult) {
    Setup();
    const auto result = matrixAdapter->CreateSparseMatrix() != nullptr;
    TearDown();
    return result == expectedResult;
}

bool MatrixAdapterTestFixture::TestCreatePhylipFile(const bool expectedResult) {
    Setup();
    const auto result = matrixAdapter->CreatePhylipFile("");
    TearDown();
    return result == expectedResult;
}

bool MatrixAdapterTestFixture::TestDistanceMatrixToSquareMatrix(const int expectedResult) {
    Setup();
    const auto result = matrixAdapter->DistanceMatrixToSquareMatrix();
    TearDown();
    return result.size() == expectedResult;
}

bool MatrixAdapterTestFixture::TestGetListVector(const bool createSparseFirst, const bool expectedResult) {
    Setup();
    if(createSparseFirst) matrixAdapter->CreateSparseMatrix();
    const auto result = matrixAdapter->GetListVector()->size() > 0;
    TearDown();
    return result == expectedResult;
}

bool MatrixAdapterTestFixture::TestGetSparseDistanceMatrix(const bool createSparseFirst, const bool expectedResult) {
    Setup();
    if(createSparseFirst) matrixAdapter->CreateSparseMatrix();
    const auto result = !matrixAdapter->GetSpareDistanceMatrix()->seqVec.empty();
    TearDown();
    return result == expectedResult;
}

MatrixAdapterTestFixture::~MatrixAdapterTestFixture() {
    delete matrixAdapter;
}

void MatrixAdapterTestFixture::Setup() {
    matrixAdapter = new MatrixAdapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25},
        0.2);
}

void MatrixAdapterTestFixture::TearDown() {
    matrixAdapter = nullptr;
}
