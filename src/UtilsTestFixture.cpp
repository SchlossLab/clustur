//
// Created by Gregory Johnson on 5/8/24.
//

#include "Tests/UtilsTestFixture.h"

bool UtilsTestFixture::TestMothurRandomShufflesRandomizesData(std::vector<int> listToRandomize,
                                                              const bool expectedResult) {
    Setup();
    //Deep copy
    const size_t count = listToRandomize.size();
    std::vector<int> copiedList(count);
    for(int i = 0; i < count; i++) {
        copiedList[i] = listToRandomize[i];
    }
    utils->mothurRandomShuffle(listToRandomize);
    TearDown();
    const bool result =  listToRandomize != copiedList;
    return result == expectedResult;
}

bool UtilsTestFixture::TestMotherRandomShuffleOverloadRandomizesData(std::vector<std::string> listToRandomize,
    const bool expectedResult) {
    Setup();
    //Deep copy
    const size_t count = listToRandomize.size();
    std::vector<std::string> copiedList(count);
    for(int i = 0; i < count; i++) {
        copiedList[i] = listToRandomize[i];
    }
    utils->mothurRandomShuffle(listToRandomize);
    TearDown();
    const bool result =  listToRandomize != copiedList;
    return result == expectedResult;
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
// Float, Double and int conversions return false if there is no error, so we invert them
bool UtilsTestFixture::TestMothurConvertStringToIntIsCorrect(const std::string &item, int number,
    const bool expectedResult) {
    Setup();
    const bool result = !utils->mothurConvert(item, number);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestMothurConvertStringToFloatIsCorrect(const std::string &item, float number,
    const bool expectedResult) {
    Setup();
    const bool result = !utils->mothurConvert(item, number);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestMothurConvertStringToDoubleIsCorrect(const std::string &item, double number,
    const bool expectedResult) {
    Setup();
    const bool result = !utils->mothurConvert(item, number);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestConvertWorksCorrectly(const std::string &str, const int expectedResult) {
    Setup();
    int num;
    utils->convert<int>(str, num);
    TearDown();
    return num == expectedResult;
}

bool UtilsTestFixture::TestGetSimpleLabelReturnsCorrectData(const std::string &label, const std::string &expectedResult) {
    Setup();
    const auto result = utils->getSimpleLabel(label);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestGetLabelTagsReturnsCorrectData(const std::string &label, const std::string &expectedResult) {
    Setup();
    const auto result = utils->getLabelTag(label);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestIsNumericIsCorrect(const std::string &stringToCheck, const bool expectedResult) {
    Setup();
    const auto result = utils->isNumeric1(stringToCheck);
    TearDown();
    return result == expectedResult;
}

bool UtilsTestFixture::TestSplitAtCommaReturnsCorrectData(std::string &stringToSplit,
    std::vector<std::string> container, const std::vector<std::string>& expectedResult) {
    Setup();
    utils->splitAtComma(stringToSplit, container);
    TearDown();
    return container == expectedResult;

}
void UtilsTestFixture::Setup() {
    utils = std::make_unique<Utils>();
}

void UtilsTestFixture::TearDown() {
    utils.reset(nullptr);
}
