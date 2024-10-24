//
// Created by gregj on 7/9/2024.
//

#include "Tests/PhylipReaderTestFixture.h"

bool PhylipReaderTestFixture::TestReadPhylipFile(const std::string &file, const int expectedResult) {
    Setup();
    reader->Read(file);
    const int result = reader->GetListVector()->getNumSeqs();
    TearDown();
    return result == expectedResult;
}

bool PhylipReaderTestFixture::TestReadPhylipFileFromRowData(const std::vector<RowData> &rowData, const int expectedResult) {
    Setup();
    reader->ReadRowDataMatrix(rowData);
    const int result = reader->GetListVector()->getNumSeqs();
    TearDown();
    return result == expectedResult;
}

bool PhylipReaderTestFixture::TestReadPhylipFileToRowData(const std::string& filePath,
    const std::vector<RowData>& expectedResult) {
    Setup();
    const auto result = reader->ReadToRowData(filePath);
    TearDown();
    return result.size() == expectedResult.size();
}

bool PhylipReaderTestFixture::TestGetDistanceMatrix(const std::vector<RowData> &rowData, const bool expectedResult) {
    Setup();
    reader->ReadRowDataMatrix(rowData);
    const auto result = !reader->GetSparseMatrix()->seqVec.empty();
    TearDown();
    return result == expectedResult;

}

bool PhylipReaderTestFixture::TestGetListVector(const std::vector<RowData> &rowData, const int expectedResult) {
    Setup();
    reader->ReadRowDataMatrix(rowData);
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
