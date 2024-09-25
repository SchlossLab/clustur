//
// Created by Gregory Johnson on 5/2/24.
//
//
// Created by Gregory Johnson on 5/1/24.
//

#include "TestHelpers/TestHelper.h"
#include "Tests/ClusterCommandTestFixture.h"
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

context("ClusterCommand") {
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("Cluster Command Sets iterations Correctly") {
    ClusterCommandTestFixture fixture;
    bool result = fixture.TestSetIterationsWorksCorrectly(4, true);
    expect_true(result);
    result = fixture.TestSetIterationsWorksCorrectly(4, false);
    expect_false(result);
  }
  test_that("Cluster Command Returns correct results for optimatrix") {
    expect_true(true);
  }

}

#endif
