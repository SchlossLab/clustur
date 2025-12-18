//
// Created by gregj on 7/9/2024.
//

#include "Tests/PhylipReaderTestFixture.h"

bool PhylipReaderTestFixture::TestReadPhylipFile(const Rcpp::DataFrame& df,
    const std::string &file, const bool expectedResult) {
    Setup();
    CountTableAdapter countTable;
    countTable.CreateDataFrameMap(df);
    reader = new ReadPhylipMatrix(countTable, 0.2, false);
    const bool result = reader->Read(file);
    TearDown();
    return result == expectedResult;
}


bool PhylipReaderTestFixture::TestGetSparseMatrix(const Rcpp::DataFrame& df,
    const std::string &file, const bool expectedResult) {
    Setup();
    CountTableAdapter countTable;
    countTable.CreateDataFrameMap(df);
    reader = new ReadPhylipMatrix(countTable, 0.2, false);
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
    CountTableAdapter countTable;
    countTable.CreateDataFrameMap(df);
    reader = new ReadPhylipMatrix(countTable, 0.2, false);
    reader->Read(file);
    const ListVector* listVector = reader->GetListVector();
    const int result = listVector->getNumSeqs();
    delete listVector;
    TearDown();
    return result == expectedResult;
}


void PhylipReaderTestFixture::Setup() {
    // Do nothing
}

void PhylipReaderTestFixture::TearDown() {
    delete reader;
}
