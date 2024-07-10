//
// Created by gregj on 7/9/2024.
//

#include "PhylipReaderTestFixture.h"

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

bool PhylipReaderTestFixture::TestGetDistanceMatix(const std::vector<RowData> &rowData, const bool expectedResult) {
    Setup();
    reader->read(rowData);
    const auto result = reader->getDMatrix()->seqVec.empty();
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
    reader = new ReadPhylipMatrix(0.2);
}

void PhylipReaderTestFixture::TearDown() {
    delete reader;
}
