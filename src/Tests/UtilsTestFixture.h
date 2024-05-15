//
// Created by Gregory Johnson on 5/8/24.
//

#ifndef UTILSTESTFIXTURE_H
#define UTILSTESTFIXTURE_H
#include <vector>
#include <memory>
#include <string>
#include "TestFixture.h"
#include "../MothurDependencies/Utils.h"


class UtilsTestFixture : TestFixture {
public:
    bool TestMothurRandomShufflesRandomizesData(std::vector<int> listToRandomize,
    bool expectedResult);
    bool TestMotherRandomShuffleOverloadRandomizesData(std::vector<std::string> listToRandomize,
        bool expectedResult);
    bool TestGetNumNamesReturnsCorrectNames(const std::string& namesToGet, int expectedResult);
    bool TestGetOTUNamesReturnsCorrectNames(std::vector<std::string> currentLabels, int numBins,
        const std::string& tagHeader, int expectedResult);
    bool TestMothurConvertStringToIntIsCorrect(const std::string& item, int number, bool expectedResult);
    bool TestMothurConvertStringToFloatIsCorrect(const std::string& item, float number, bool expectedResult);
    bool TestMothurConvertStringToDoubleIsCorrect(const std::string& item, double number, bool expectedResult);
    bool TestConvertWorksCorrectly(const std::string& str, int expectedResult);
    bool TestGetSimpleLabelReturnsCorrectData(const std::string& label, const std::string& expectedResult);
    bool TestGetLabelTagsReturnsCorrectData(const std::string& label, const std::string& expectedResult);
    bool TestIsNumericIsCorrect(const std::string& stringToCheck, bool expectedResult);
    bool TestSplitAtCommaReturnsCorrectData(std::string& stringToSplit, std::vector<std::string> container,
        const std::vector<std::string>& expectedResult);
    bool TestGetOTUNamesReturnsOTUNames(const std::vector<std::string>& currentLabels, int numBins,
        const std::string& tagHeader, int expectedResult);
private:
    std::unique_ptr<Utils> utils;
    void Setup() override;
    void TearDown() override;
};



#endif //UTILSTESTFIXTURE_H
