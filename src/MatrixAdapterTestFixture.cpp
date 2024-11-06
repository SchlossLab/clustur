//
// Created by Gregory Johnson on 7/12/24.
//

#include "Tests/MatrixAdapterTestFixture.h"

bool MatrixAdapterTestFixture::TestCreateSparseMatrix(const size_t expectedResult) {
    Setup();
    const auto result = matrixAdapter->CreateSparseMatrix().seqVec.size();
    TearDown();
    return result == expectedResult;
}

bool MatrixAdapterTestFixture::TestCreatePhylipFile(const bool expectedResult) {
    Setup();
    const auto result = matrixAdapter->CreatePhylipFile("");
    TearDown();
    return result == expectedResult;
}

bool MatrixAdapterTestFixture::TestGetListVector(const bool createSparseFirst, const bool expectedResult) {
    Setup();
    if(createSparseFirst) matrixAdapter->CreateSparseMatrix();
    const auto result = matrixAdapter->CreateListVector().size() > 0;
    TearDown();
    return result == expectedResult;
}

MatrixAdapterTestFixture::~MatrixAdapterTestFixture() {
    delete matrixAdapter;
}

void MatrixAdapterTestFixture::Setup() {
    const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
    const std::vector<double> total{10, 20, 30, 40, 50, 60};
    const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
        Rcpp::Named("Representative Sequence") = compounds,
        Rcpp::Named("total") = total,
        Rcpp::Named("nogroup") = total);
    CountTableAdapter countTable;
    countTable.CreateDataFrameMap(dataframe);
    matrixAdapter = new MatrixAdapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25},
        0.2, false, countTable);
}

void MatrixAdapterTestFixture::TearDown() {
    matrixAdapter = nullptr;
}
