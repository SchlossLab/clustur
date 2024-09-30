//
// Created by gregj on 7/9/2024.
//

#include "Tests/PhylipReaderTestFixture.h"

bool PhylipReaderTestFixture::TestReadPhylipFile(const std::string &file, const bool expectedResult) {
    Setup();
    const auto result = reader->read(file);
    TearDown();
    return result == expectedResult;
}

bool PhylipReaderTestFixture::TestReadPhylipFile(const std::vector<RowData> &rowData, const bool expectedResult) {
    Setup();
    const auto result = reader->read(rowData);
    TearDown();
    return result == expectedResult;
}

bool PhylipReaderTestFixture::TestReadPhylipFileToRowData(const std::string& filePath,
    const std::vector<RowData>& expectedResult) {
    Setup();
    const auto result = reader->readToRowData(filePath);
    TearDown();
    return !result.empty() && !expectedResult.empty();
}

bool PhylipReaderTestFixture::TestGetDistanceMatrix(const std::vector<RowData> &rowData, const bool expectedResult) {
    Setup();
    reader->read(rowData);
    const auto result = !reader->getDMatrix()->seqVec.empty();
    TearDown();
    return result == expectedResult;

}

bool PhylipReaderTestFixture::TestGetListVector(const std::vector<RowData> &rowData, const bool expectedResult) {
    Setup();
    reader->read(rowData);
    const auto result = reader->getListVector()->size() > 0;
    TearDown();
    return result == expectedResult;
}


void PhylipReaderTestFixture::Setup() {
    reader = new ReadPhylipMatrix(0.2, false);
}

void PhylipReaderTestFixture::TearDown() {
    delete reader;
}
