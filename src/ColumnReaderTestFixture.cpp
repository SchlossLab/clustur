//
// Created by Gregory Johnson on 9/16/24.
//

#include "Tests/ColumnReaderTestFixture.h"

bool ColumnReaderTestFixture::TestReadColumnFile(const std::string &file, const Rcpp::DataFrame& df,
                                                 const int expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(file);
    const int result = reader->GetListVector()->getNumSeqs();
    TearDown();
    return result == expectedResult;
}

bool ColumnReaderTestFixture::TestGetDistanceMatrix(const std::string& filePath,
    const Rcpp::DataFrame& df, const bool expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(filePath);
    const auto result = !reader->GetSparseMatrix()->seqVec.empty();
    TearDown();
    return result == expectedResult;
}

bool ColumnReaderTestFixture::TestGetListVector(const std::string& filePath,
    const Rcpp::DataFrame& df, const int expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(filePath);
    const auto result = reader->GetListVector()->getNumSeqs();
    TearDown();
    return result == expectedResult;
}

void ColumnReaderTestFixture::Setup() {
    reader = new ColumnDistanceMatrixReader(0.2, false);
}

void ColumnReaderTestFixture::TearDown() {
    delete reader;
}
