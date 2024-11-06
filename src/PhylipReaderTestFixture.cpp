//
// Created by gregj on 7/9/2024.
//

#include "Tests/PhylipReaderTestFixture.h"

bool PhylipReaderTestFixture::TestReadPhylipFile(const Rcpp::DataFrame& df,
    const std::string &file, const bool expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    const bool result = reader->Read(file);
    TearDown();
    return result == expectedResult;
}


bool PhylipReaderTestFixture::TestGetSparseMatrix(const Rcpp::DataFrame& df,
    const std::string &file, const bool expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(file);
    const auto result = !reader->GetSparseMatrix()->seqVec.empty();
    TearDown();
    return result == expectedResult;

}

bool PhylipReaderTestFixture::TestGetListVector(const Rcpp::DataFrame& df,
    const std::string &file, const int expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(file);
    const int result = reader->GetListVector()->getNumSeqs();
    TearDown();
    return result == expectedResult;
}


void PhylipReaderTestFixture::Setup() {
    reader = new ReadPhylipMatrix(0.2, false);
}

void PhylipReaderTestFixture::TearDown() {
    delete reader;
}
