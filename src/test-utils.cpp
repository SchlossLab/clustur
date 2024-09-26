//
// Created by Gregory Johnson on 5/14/24.
//
#include "TestHelpers/TestHelper.h"
#include "Tests/UtilsTestFixture.h"
#if DEBUG_RCPP
#include <testthat.h>
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
//  OptiMatrix* ConvertToOptimatrix(const std::vector<int>
//        &xPosition, const std::vector<int>& yPosition, const std::vector<double>& data, int rowSize, int colSize);
// void Check()
// {
//   OptimatrixAdapterTestFixture fixture;
//   std::vector<int> xVals = std::vector<int>{1,1,1,2,2,3,5};
//   std::vector<int> yVals = std::vector<int>{2,3,5,3,5,5,5};
//   std::vector<double> data = std::vector<double>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
//   bool result = fixture.TestOptimatrixClosenessReturnsCorrectValue(xVals, yVals, data, 4);
// }

context("Utils Testing") {
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("Utils Converts properly") {
    UtilsTestFixture testFixture;
    const std::string number = "123";
    bool result = testFixture.TestConvertWorksCorrectly(number, 123);
    expect_true(result);
  }
  test_that("Utils is numeric works properly") {
    UtilsTestFixture testFixture;
    const std::string number = "123";
    bool result = testFixture.TestIsNumericIsCorrect(number, true);
    expect_true(result);
  }
  test_that("Utils MothurRandomShuffle works as expected") {
    UtilsTestFixture testFixture;
    const std::vector<int> exampleIntConainer = {1,2,3,4};
    bool result = testFixture.TestMothurRandomShufflesRandomizesData(exampleIntConainer, true);
    expect_true(result);
    result = testFixture.TestMothurRandomShufflesRandomizesData({}, false);
    expect_true(result);
  }
  test_that("Utils Get Labels returns labels correctly") {
    UtilsTestFixture testFixture;
    bool result = testFixture.TestGetLabelTagsReturnsCorrectData("OTU300", "OTU");
    expect_true(result);
  }
  test_that("Utils Get Simple Labels returns the number correctly") {
    UtilsTestFixture testFixture;
    bool result = testFixture.TestGetSimpleLabelReturnsCorrectData("OTU200", "200");
    expect_true(result);
  }
  test_that("Utils shuffle, shuffles string list correctly") {
    UtilsTestFixture testFixture;
    const std::vector<std::string> exampleContainer = {"1", "2", "3", "4"};
    bool result = testFixture.TestMotherRandomShuffleOverloadRandomizesData(exampleContainer, true);
    expect_true(result);
  }
  test_that("Utils splits at commas, returns a vector of all the data split") {
    UtilsTestFixture testFixture;
    const std::vector<std::string> exampleContainer = {"1", "2", "3", "4"};
    std::string exampleData = "1,2,3,4";
    bool result = testFixture.TestSplitAtCommaReturnsCorrectData(exampleData, {}, exampleContainer);
    expect_true(result);
  }
  test_that("Utils Get Number of names returns the correct number") {
    UtilsTestFixture testFixture;
    bool result = testFixture.TestGetNumNamesReturnsCorrectNames("1,2,4", 3);
    expect_true(result);
  }
  test_that("Utils Converts Doubles properly") {
    UtilsTestFixture testFixture;
    bool result = testFixture.TestMothurConvertStringToDoubleIsCorrect("2.00", 2, true);
    expect_true(result);
    result = testFixture.TestMothurConvertStringToDoubleIsCorrect("==", 2, false);
    expect_true(result);

  }
  test_that("Utils Converts Floats Properly") {
    UtilsTestFixture testFixture;
    bool result = testFixture.TestMothurConvertStringToFloatIsCorrect("2.00", 2, true);
    expect_true(result);
    result = testFixture.TestMothurConvertStringToFloatIsCorrect("==", 2, false);
    expect_true(result);
  }
  test_that("Utils Converts Ints Properly") {
    UtilsTestFixture testFixture;
    bool result = testFixture.TestMothurConvertStringToIntIsCorrect("2", 2, true);
    expect_true(result);
    result = testFixture.TestMothurConvertStringToIntIsCorrect("==", 2, false);
    expect_true(result);
  }
  test_that("Get OTU Names returns the correct number of OTU Names") {
    UtilsTestFixture testFixture;
    bool result = testFixture.TestGetOTUNamesReturnsCorrectNames({"OTU12", "OTU13"}, 4, "OTU", 4);
    expect_true(result);
    result = testFixture.TestGetOTUNamesReturnsCorrectNames({"OTU12", "OTU13"}, 2, "OTU", 0);
    expect_true(result);
  }
}

#endif
