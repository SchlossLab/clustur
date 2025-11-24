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
    const SparseDistanceMatrix* sparseMatrix = reader->GetSparseMatrix();
    const auto result = !sparseMatrix->seqVec.empty();
    delete sparseMatrix;
    TearDown();
    return result == expectedResult;

}

bool PhylipReaderTestFixture::TestGetListVector(const Rcpp::DataFrame& df,
    const std::string &file, const int expectedResult) {
    Setup();
    reader->CreateCountTableAdapter(df);
    reader->Read(file);
    const ListVector* listVector = reader->GetListVector();
    const int result = listVector->getNumSeqs();
    delete listVector;
    TearDown();
    return result == expectedResult;
}


void PhylipReaderTestFixture::Setup() {
    reader = new ReadPhylipMatrix(0.2, false);
}

void PhylipReaderTestFixture::TearDown() {
    delete reader;
}
