//
// Created by Gregory Johnson on 5/7/24.
//

#ifndef OPTIDATATESTFIXTURE_H
#define OPTIDATATESTFIXTURE_H
#include <set>
#include <string>
#include <memory>
#include "TestFixture.h"
#include "../MothurDependencies/ListVector.h"
#include "../MothurDependencies/OptiData.h"


class OptiDataTestFixture : TestFixture {
public:
    bool TestGetCloseSequencesReturnCorrectData(long long index, const std::unordered_set<long long>& expectedResults);
    bool TestIsCloseReturnsIfDataIsClose(long long index, long long indexToFind, bool expectedResult);
    bool TestGetNumCloseReturnsCorrectData(long long index, long long expectedResult);
    bool TestGetNameReturnsTheProperName(long long index, const std::string& expectedResult);
    bool TestGetNumSeqsReturnsCorrectData(long long expectedResult);
    bool TestGetNumSingletonsReturnsCorrectData(long long expectedResult);
    bool TestGetListSingleReturnsCorrectData(const int expectedResult);
    bool TestGetCloseRefSeqsReturnsCorrectData(long long index, const std::unordered_set<long long>& expectedResult);
    bool TestIsCloseFitReturnsCorrectData(long long indexOne, long long indexTwo, bool expected, bool expectedResult);

private:
    OptiData* optiData = nullptr;
    void Setup() override;
    void TearDown() override;
};



#endif //OPTIDATATESTFIXTURE_H
