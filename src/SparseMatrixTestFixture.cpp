//
// Created by gregj on 7/9/2024.
//

#include "Tests/SparseMatrixTestFixture.h"

#include "Adapters/MatrixAdapter.h"

bool SparseMatrixTestFixture::TestGetSmallestCell(unsigned long index, const float expectedResult) {
    Setup();
    const unsigned long smallestCell = sparseDistanceMatrix->getSmallestCell(index);
    float result = -1000;
    for(const auto &seq:  sparseDistanceMatrix->seqVec[index]) {
        if(seq.index == smallestCell) {
            result = seq.dist;
            break;
        }
    }
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
    unsigned long currentSize = 0;
    if(row < sparseDistanceMatrix->seqVec.max_size())
       currentSize = sparseDistanceMatrix->seqVec[row].size();
    sparseDistanceMatrix->rmCell(row, col);
    const bool result = sparseDistanceMatrix->seqVec[row].size() < currentSize;
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestUpdateCellCompliment(const unsigned long row, const unsigned long col,
    const bool expectedResult) {
    Setup();
    sparseDistanceMatrix->updateCellCompliment(row, col);
    const bool result =  sparseDistanceMatrix->seqVec[3][2].dist == sparseDistanceMatrix->seqVec[row][col].dist;
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestResize(const unsigned long size, const long expectedResult) {
    Setup();
    const size_t currentSize = sparseDistanceMatrix->seqVec.size();
    sparseDistanceMatrix->resize(size);
    const unsigned long result = sparseDistanceMatrix->seqVec.size();
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

bool SparseMatrixTestFixture::TestAddCell(const unsigned long row, const PDistCell& cell, const bool expectedResult) {
    Setup();
    const size_t size = sparseDistanceMatrix->seqVec[row].size();
    sparseDistanceMatrix->addCell(row, cell);
    const bool result = sparseDistanceMatrix->seqVec[row].size() > size;
    TearDown();
    return result == expectedResult;

}

bool SparseMatrixTestFixture::TestAddCellSorted(const unsigned long row, const PDistCell& cell, const bool expectedResult) {
    Setup();
    sparseDistanceMatrix->addCellSorted(row, cell);
    const size_t size = sparseDistanceMatrix->seqVec[row].size() - 1;
    const auto currentCell = sparseDistanceMatrix->seqVec[row][size - cell.index];
    const bool result = currentCell.dist == cell.dist && currentCell.index == cell.index;
    TearDown();
    return result == expectedResult;
}

bool SparseMatrixTestFixture::TestPrint(const bool clear, const bool expectedResult) {
    Setup();
    if(clear)
        sparseDistanceMatrix->clear();
    const bool result = sparseDistanceMatrix->print();
    TearDown();
    return result == expectedResult;
}


void SparseMatrixTestFixture::Setup() {
    MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false);
    sparseDistanceMatrix = adapter.CreateSparseMatrix();
}

void SparseMatrixTestFixture::TearDown() {
    delete sparseDistanceMatrix;
}
