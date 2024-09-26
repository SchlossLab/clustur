//
// Created by Gregory Johnson on 5/1/24.
//

#include "TestHelpers/TestHelper.h"
#if DEBUG_RCPP
#include <testthat.h>
#include "Tests/ListVectorTestFixture.h"
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
context("ListVectorTest") {
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("List Vector Returns proper amount of bins") {
    std::vector<std::string> mockData {"1", "2", "3", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorReturnsCorrectNumberOfBins(mockData, 10);
    expect_true(result);
    // Test for the negative
    result = testFixture.TestListVectorReturnsCorrectNumberOfBins(mockData, 9);
    expect_false(result);
  }
  test_that("List Vector Returns the correct amount of numSeqs") {
    std::vector<std::string> mockData {"1", "2", "3", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorReturnsCorrectNumSeqs(mockData, 10);
    mockData = {"1", "2", "3", "4", "5", "6", "7,2", "8,4","9", "10"};
    result = testFixture.TestListVectorReturnsCorrectNumSeqs(mockData, 12);
    expect_true(result);
    // Test the negative
    result = testFixture.TestListVectorReturnsCorrectNumSeqs(mockData, 10);
    expect_false(result);
  }
  test_that("List Vector Returns the correct MaxRank") {
    std::vector<std::string> mockData {"1", "2,2", "3,3,4", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorReturnsCorrectMaxRank(mockData, 3);
    expect_true(result);
    mockData = {"1", "2,2", "3,3", "4", "5", "6", "7", "8","9", "10"};
    result = testFixture.TestListVectorReturnsCorrectMaxRank(mockData, 2);
    expect_true(result);
    result = testFixture.TestListVectorReturnsCorrectMaxRank(mockData, 1);
    expect_false(result);
  }
  test_that("List Vector Returns the correct GetValue Based on index") {
    std::vector<std::string> mockData {"1", "2,2", "3,3,4", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorReturnsCorrectGetValue(mockData, 0, "1");
    expect_true(result);
    result = testFixture.TestListVectorReturnsCorrectGetValue(mockData, 1, "2,2");
    expect_true(result);
    result = testFixture.TestListVectorReturnsCorrectGetValue(mockData, -1, "");
    expect_true(result);
    result = testFixture.TestListVectorReturnsCorrectGetValue(mockData, 3, "");;
    expect_false(result);
  }
  test_that("List Vector Returns the correct label values") {
    std::vector<std::string> mockData {"1", "2,2", "3,3,4", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorReturnsCorrectGetLabelsValue(mockData, 10);
    expect_true(result);
  }
  test_that("List Vector Returns correct OTU size") {
    std::vector<std::string> mockData {"1", "2,2", "3,3,4", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorReturnsCorrectGetOtuNamesSize(mockData, 2, 2);
    expect_true(result);
  }
  test_that("List Vector Returns OTU Names") {
    std::vector<std::string> mockData {"1", "2,2", "3,3,4", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorReturnsCorrectGetOtuNames(mockData, 1, "02");
    expect_true(result);
  }
  test_that("List Vector Returns the data given a index") {
    std::vector<std::string> mockData {"1", "2,2", "3,3,4", "4", "5", "6", "7", "8","9", "10"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorPushBackSetsDataCorrectly(mockData, 1, "2,2");
    expect_true(result);
  }
  test_that("List Vector Prints the proper data") {
    std::vector<std::string> mockData {"1", "2,2", "3"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorPrintDisplaysDataCorrectly(mockData, "\tOtu1\tOtu2\tOtu3\n\t2,2\t1\t3\n");
    expect_true(result);
  }
  test_that("List Vector Properly Sets Its Labels") {
    std::vector<std::string> mockData {"1", "2,2", "3"};
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorSetsLabelsCorrectly(mockData, 3);
    expect_true(result);
  }
  test_that("Set PrintedLabel Functions Correctly") {
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestListVectorSetsPrintedLabelsCorrectly(true, true);
    expect_true(result);
    result = testFixture.TestListVectorSetsPrintedLabelsCorrectly(false, true);
    expect_false(result);
  }
  test_that("Create Dataframe from list works") {
    ListVectorTestFixture testFixture;
    bool result = testFixture.TestCreateDataFrameFromList("0.0001", true);
    expect_true(result);
    result = testFixture.TestCreateDataFrameFromList("0.0001", false);
    expect_false(result);
  }



}

#endif
