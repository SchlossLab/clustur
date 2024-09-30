//
// Created by Gregory Johnson on 9/16/24.
//

#include "Tests/ColumnReaderTestFixture.h"

bool ColumnReaderTestFixture::TestReadColumnFile(const std::string &file, const CountTableAdapter& adapter,
                                                 const bool expectedResult) {
    Setup();
    const auto result = reader->Read(adapter, file);
    TearDown();
    return result == expectedResult;
}


bool ColumnReaderTestFixture::TestReadColumnFileToRowData(const std::string &filePath, const CountTableAdapter& adapter,
    const std::vector<RowData> &expectedResult) {
    Setup();
    const auto result = reader->readToRowData(adapter, filePath);
    TearDown();
    return !result.empty() && !expectedResult.empty();
}

bool ColumnReaderTestFixture::TestGetDistanceMatrix(const std::string& filePath,
    const CountTableAdapter& adapter, const bool expectedResult) {
    Setup();
    reader->Read(adapter, filePath);
    const auto result = !reader->GetSparseDataMatrix()->seqVec.empty();
    TearDown();
    return result == expectedResult;
}

bool ColumnReaderTestFixture::TestGetListVector(const std::string& filePath,
     const CountTableAdapter& adapter, const bool expectedResult) {
    Setup();
    reader->Read(adapter, filePath);
    const auto result = reader->GetListVector()->size() > 0;
    TearDown();
    return result == expectedResult;
}

void ColumnReaderTestFixture::Setup() {
    reader = new ColumnDistanceMatrixReader(0.2, false);
}

void ColumnReaderTestFixture::TearDown() {
    delete reader;
}
