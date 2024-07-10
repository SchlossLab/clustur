//
// Created by gregj on 7/9/2024.
//

#include "SparseMatrixTestFixture.h"

#include "../Adapters/MatrixAdapter.h"

bool SparseMatrixTestFixture::TestGetSmallestCell(unsigned long &index, const unsigned long expectedResult) {
    Setup();
    const auto result = sparseDistanceMatrix->getSmallestCell(index);
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestGetSmallestDist(const float expectedResult) {
    Setup();
    const float result = sparseDistanceMatrix->getSmallDist();
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestHeapComparator(const PDistCell &a, const PDistCell &b,
    const bool expectedResult) {
    Setup();
    const auto result = SparseDistanceMatrix::heapComparator(a, b);
    TearDown();
    return result == expectedResult;

}

bool SparseMatrixTestFixture::TestRemoveCell(const unsigned long row, const unsigned long col,
    const bool expectedResult) {
    Setup();
    const bool result = sparseDistanceMatrix->rmCell(row, col);
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestUpdateCellCompliment(const unsigned long row, const unsigned long col,
    const bool expectedResult) {
    Setup();
    const bool result = sparseDistanceMatrix->updateCellCompliment(row, col);
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestResize(const unsigned long size, const long expectedResult) {
    Setup();
    const size_t currentSize = sparseDistanceMatrix->seqVec.size();
    sparseDistanceMatrix->resize(size);
    const bool result = sparseDistanceMatrix->seqVec.size() != currentSize;
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestClear(const bool expectedResult) {
    Setup();
    sparseDistanceMatrix->clear();
    const bool result = sparseDistanceMatrix->seqVec.empty();
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestAddCell(const unsigned long row , const PDistCell& cell, const bool expectedResult) {
    Setup();
    const size_t size = sparseDistanceMatrix->seqVec.size();
    sparseDistanceMatrix->addCell(row, cell);
    const bool result = sparseDistanceMatrix->seqVec.size() > size;
    TearDown();
    return result == expectedResult;

}

bool SparseMatrixTestFixture::TestAddCellSorted(const unsigned long row, const PDistCell& cell, const bool expectedResult) {
    Setup();
    sparseDistanceMatrix->addCellSorted(row, cell);
    const auto currentCell = sparseDistanceMatrix->seqVec[row][cell.index - 1];
    const bool result = currentCell.dist == cell.dist && currentCell.index == cell.index;
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestPrint(const bool expectedResult) {
    Setup();
    const bool result = sparseDistanceMatrix->print();
    TearDown();
    return result == expectedResult;
}

SparseMatrixTestFixture::~SparseMatrixTestFixture() {
    delete sparseDistanceMatrix;
}

void SparseMatrixTestFixture::Setup() {
    MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2);
    sparseDistanceMatrix = adapter.CreateSparseMatrix();
}

void SparseMatrixTestFixture::TearDown() {
    delete sparseDistanceMatrix;
}
