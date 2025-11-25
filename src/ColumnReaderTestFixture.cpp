//
// Created by Gregory Johnson on 9/16/24.
//

#include "Tests/ColumnReaderTestFixture.h"

bool ColumnReaderTestFixture::TestReadColumnFile(const std::string &file, const Rcpp::DataFrame& df,
                                                 const int expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(file);
    const ListVector* listVector = reader->GetListVector();
    const int result = listVector->getNumSeqs();
    delete listVector;
    TearDown();
    return result == expectedResult;
}

bool ColumnReaderTestFixture::TestGetDistanceMatrix(const std::string& filePath,
    const Rcpp::DataFrame& df, const bool expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(filePath);
    const SparseDistanceMatrix* sparseMatrix = reader->GetSparseMatrix();
    const auto result = !sparseMatrix->seqVec.empty();
    delete sparseMatrix;
    TearDown();
    return result == expectedResult;
}

bool ColumnReaderTestFixture::TestGetListVector(const std::string& filePath,
    const Rcpp::DataFrame& df, const int expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(filePath);
    const ListVector* listVector = reader->GetListVector();
    const auto result = listVector->getNumSeqs();
    delete listVector;
    TearDown();
    return result == expectedResult;
}

void ColumnReaderTestFixture::Setup() {
    reader = new ColumnDistanceMatrixReader(0.2, false);
}

void ColumnReaderTestFixture::TearDown() {
    delete reader;
}
