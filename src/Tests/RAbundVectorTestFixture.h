//
// Created by Gregory Johnson on 7/12/24.
//

#ifndef RABUNDVECTORTESTFIXTURE_H
#define RABUNDVECTORTESTFIXTURE_H
#include "TestFixture.h"
#include <string>
#include <map>
#include <vector>
#include "../MothurDependencies/RAbundVector.h"




class RAbundVectorTestFixture final : public TestFixture {
public:
    bool TestRemove(int bin , int expectedResult);
    bool TestSet(int binNumber, int newBinSize, bool expectedResult);
    bool TestGet(int index, int expectedResult);
    bool TestGetOverload(const std::vector<int> &expectedResult);
    bool TestPushBack(int binSize, bool expectedResult);
    bool TestResize(int size, bool expectedResult);
    bool TestSize(int expectedResult);
    bool TestQuicksort(bool expectedResult);
    bool TestClear(bool expectedResult);
    bool TestGetNumBins(int expectedResult);
    bool TestGetNumSeqs(int expectedResult);
    bool TestGetMaxRank(int expectedResult);
private:
    RAbundVector radbund;
    void Setup() override;
    void TearDown() override;
};



#endif //RABUNDVECTORTESTFIXTURE_H
