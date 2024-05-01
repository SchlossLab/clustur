//
// Created by Gregory Johnson on 5/1/24.
//

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
context("ListVectorTest") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("List Vector Returns proper amount of bins") {
    expect_true(true);
  }

}

#endif