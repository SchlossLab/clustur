//
// Created by gregj on 7/9/2024.
//

#ifndef SPARSEMATRIXTESTFIXTURE_H
#define SPARSEMATRIXTESTFIXTURE_H
#include "TestFixture.h"
#include <vector>
#include "../MothurDependencies/sparsedistancematrix.h"


class SparseMatrixTestFixture final : public TestFixture {
public:
    bool TestGetSmallestCell(unsigned long index, float expectedResult);
    bool TestGetSmallestDist(float expectedResult);
    bool TestHeapComparator(const PDistCell& a, const PDistCell& b, bool expectedResult);
    bool TestRemoveCell(unsigned long, unsigned long, bool expectedResult);
    bool TestUpdateCellCompliment(unsigned long, unsigned long, bool expectedResult);
    bool TestResize(unsigned long size, long expectedResult);
    bool TestClear(bool expectedResult);
    bool TestAddCell(unsigned long, const PDistCell&, bool expectedResult);
    bool TestAddCellSorted(unsigned long, const PDistCell&, bool expectedResult);
    bool TestPrint(bool clear, bool expectedResult);
private:
    void Setup() override;
    void TearDown() override;
    SparseDistanceMatrix* sparseDistanceMatrix = nullptr;
};



#endif //SPARSEMATRIXTESTFIXTURE_H
