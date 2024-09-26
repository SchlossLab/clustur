
#include "TestHelpers/TestHelper.h"
#if DEBUG_RCPP
#include <testthat.h>
#include "Tests/OptimatrixAdapterTestFixture.h"
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
context("Optimatrix Adapter Test") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("OptimatrixAdapter returns a optimatrix") {
    OptimatrixAdapterTestFixture fixture;
    bool result = fixture.TestOptimatrixReturnsCorrectValue(std::vector<int>(), std::vector<int>(), std::vector<double>());
    expect_true(result);
  }
  test_that("OptimatrixAdapter returns the correct number of Closeness Values"){
    OptimatrixAdapterTestFixture fixture;
    const auto xVals = std::vector<int>{0,0,0,1,1,2,3};
    const auto yVals = std::vector<int>{1,2,4,2,4,4,4};
    const auto data = std::vector<double>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
    bool result = fixture.TestOptimatrixClosenessReturnsCorrectValue(xVals, yVals, data, 4);
    expect_true(result);
  }
  test_that("OptimatrixAdapter returns the correct number of Singletons"){
    OptimatrixAdapterTestFixture fixture;
    const auto xVals = std::vector<int>{0,0,0,1,1,2,3};
    const auto yVals = std::vector<int>{1,2,4,2,4,4,4};
    const auto data = std::vector<double>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
    bool result = fixture.TestOptimatrixSingletonReturnsCorrectValue(xVals, yVals, data, 1);
    expect_true(result);
  }
  test_that("OptimatrixAdapter returns the correct number of Names"){
    OptimatrixAdapterTestFixture fixture;
    const auto xVals = std::vector<int>{0,0,0,1,1,2,3};
    const auto yVals = std::vector<int>{1,2,4,2,4,4,4};
    const auto data = std::vector<double>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
    bool result = fixture.TestOptimatrixNameListReturnsCorrectValue(xVals, yVals, data, 4);
    expect_true(result);
  }
  test_that("OptimatrixAdapter Convert properly converts rowdata into a OptimatrixAdapter"){
    OptimatrixAdapterTestFixture test;
    std::vector<RowData> data(3);
    data[0].name = "0";
    data[1].name = "1";
    data[2].name = "2";
    data[0].rowValues = {0,0.3,0};
    data[1].rowValues = {0.3,0,0};
    data[2].rowValues = {0.3,0,0};
    std::vector<std::set<long long>> close(2);
    close[0].insert(0);
    close[0].insert(1);
    close[1].insert(0);
    close[1].insert(1);
    bool result = test.ConvertToOptimatrixWithRowData(data, false, close);
    expect_true(result);
    result = test.ConvertToOptimatrixWithRowData(data, false, {});
    expect_false(result);
  }

}

#endif
