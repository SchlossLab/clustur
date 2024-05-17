//
// Created by Gregory Johnson on 5/1/24.
//

#ifndef LISTVECTORTESTFIXTURE_H
#define LISTVECTORTESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/ListVector.h"
#include <memory>

class ListVectorTestFixture : TestFixture {
public:
    bool TestListVectorReturnsCorrectNumberOfBins(const std::vector<std::string>& mockListOfSequences,
        int expectedResult);
    bool TestListVectorReturnsCorrectNumSeqs(const std::vector<std::string>& mockListOfSequences,
        int expectedResult);
    bool TestListVectorReturnsCorrectMaxRank(const std::vector<std::string>& mockListOfSequences,
        int expectedResult);
    bool TestListVectorReturnsCorrectGetValue(const std::vector<std::string>& mockListOfSequences,
        int index, const std::string& expectedResult);
    bool TestListVectorReturnsCorrectGetLabelsValue(const std::vector<std::string>& mockListOfSequences,
        int expectedResult);
    bool TestListVectorReturnsCorrectGetOtuNamesSize(const std::vector<std::string>& mockListOfSequences,
        int binToTest, int expectedResult);
    bool TestListVectorReturnsCorrectGetOtuNames(const std::vector<std::string>& mockListOfSequences,
        int binToTest, const std::string& expectedResult);
    bool TestListVectorSetsLabelsCorrectly(const std::vector<std::string>& mockListOfLabels,
        int expectedResult);
    bool TestListVectorSetsPrintedLabelsCorrectly(bool mockValue, bool expectedResult);
    bool TestListVectorPushBackSetsDataCorrectly(const std::vector<std::string>& mockListOfSequences,
        int index, const std::string& expectedResult);
    bool TestListVectorPrintDisplaysDataCorrectly(const std::vector<std::string>& mockListOfSequences,
        const std::string& expectedResult);
private:
    void Setup() override;
    void TearDown() override;
    std::unique_ptr<ListVector> listVector;
};



#endif //LISTVECTORTESTFIXTURE_H
