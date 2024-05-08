//
// Created by Gregory Johnson on 5/8/24.
//

#include "UtilsTestFixture.h"

bool UtilsTestFixture::TestMothurRandomShufflesRandomizesData(std::vector<int> listToRandomize,
    const std::vector<int> &expectedResult) {
    Setup();
    utils->mothurRandomShuffle(listToRandomize);
    TearDown();
    return listToRandomize == expectedResult;
}

bool UtilsTestFixture::TestMotherRandomShuffleOverloadRandomizesData(std::vector<std::string> listToRandomize,
    const std::vector<std::string> &expectedResult) {
    Setup();
    utils->mothurRandomShuffle(listToRandomize);
    TearDown();
    return listToRandomize == expectedResult;
}

bool UtilsTestFixture::TestGetNumNamesReturnsCorrectNames(const std::string& namesToGet, const int expectedResult) {
    Setup();
    const int result = utils->getNumNames(namesToGet);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestGetOTUNamesReturnsCorrectNames(std::vector<std::string> currentLabels, const int numBins,
    const std::string &tagHeader, const int expectedResult) {
    Setup();
    const int result = utils->getOTUNames(currentLabels, numBins, tagHeader);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestMothurConvertStringToIntIsCorrect(const std::string &item, int number, bool expectedResult) {
    Setup();
    const bool result = utils->mothurConvert(item, number);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestMothurConvertStringToFloatIsCorrect(const std::string &item, const float number,
    bool expectedResult) {
}

bool UtilsTestFixture::TestMothurConvertStringToDoubleIsCorrect(const std::string &item, const double number,
    bool expectedResult) {
}

bool UtilsTestFixture::TestConvertWorksCorrectly(const std::string &str, bool failIfLeftoverChars) {
}

bool UtilsTestFixture::
TestGetSimpleLabelReturnsCorrectData(const std::string &label, const std::string &expectedResult) {
}

bool UtilsTestFixture::TestGetLabelTagsReturnsCorrectData(const std::string &label, const std::string &expectedResult) {
}

bool UtilsTestFixture::TestIsNumericIsCorrect(const std::string &stringToCheck, bool expectedResult) {
}

bool UtilsTestFixture::TestSplitAtCommaReturnsCorrectData(const std::string &stringToSplit,
    std::vector<std::string> container) {
}

bool UtilsTestFixture::TestSplitAtCommaOverloadReturnsCorrectData(std::string prefix, std::string suffix,
    const std::string &expectedResultPrefix, const std::string &expectedResultSuffix) {
}

bool UtilsTestFixture::TestSplitAtCommaIntOVerloadReturnsCorrectData(std::string stringToSplit,
    std::vector<int> container, std::vector<int> containerResultExpected, std::string stringToSplitExpectedResult) {
}

void UtilsTestFixture::Setup() {
    utils = std::make_unique<Utils>();
}

void UtilsTestFixture::TearDown() {
}
