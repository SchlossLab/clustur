//
// Created by Gregory Johnson on 5/6/24.
//
//
// Created by Gregory Johnson on 5/1/24.
//

#include "TestHelpers/TestHelper.h"
#include "Tests/OptiClusterTestFixture.h"
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
// }
context("Opticluster test") {
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("Opticluster initializes properly") {
    OptiClusterTestFixture testFixture;
    bool res =  testFixture.OptiClusterInitalizesProperly(1, false, "", 0);
    expect_true(res);
    res =  testFixture.OptiClusterInitalizesProperly(1, false, "", 1);
    expect_false(res);

  }
  test_that("Opticluster Updates properly") {
    OptiClusterTestFixture testFixture;
    bool res = testFixture.OptiClusterUpdatesCorrectly(0.5, false);
    expect_true(res);
    res = testFixture.OptiClusterUpdatesCorrectly(0.5, true);
    expect_false(res);

  }
  test_that("Opticluster returns proper stats") {
    OptiClusterTestFixture testFixture;
    bool result = testFixture.OpticlusterReturnsTheCorrectStats(1,1,1,1, {0,0,0,0,1,0,0,0});
    expect_true(result);
    result = testFixture.OpticlusterReturnsTheCorrectStats(1,1,1,1, {0,0,0,0,1,1,0,0});
    expect_false(result);
  }
  test_that("Opticluster Gets CloseFarCounts properly") {
    OptiClusterTestFixture testFixture;
    bool res = testFixture.OptiClusterGetsTheCorrectCloseFarCounts(1, 2, {0,0});
    expect_true(res);
    res = testFixture.OptiClusterGetsTheCorrectCloseFarCounts(1, 2, {0,1});
    expect_false(res);
  }
  test_that("Opticluster Returns ListVector properly") {
    OptiClusterTestFixture testFixture;
    bool res = testFixture.OptiClusterReturnsTheCorrectList(1);
    expect_true(res);
    res = testFixture.OptiClusterReturnsTheCorrectList(0);
    expect_false(res);
  }
  test_that("Opticluster Gets the correct tags") {
    OptiClusterTestFixture testFixture;
    bool res = testFixture.OptiClusterGetsTheCorrectTag("opti_mcc");
    expect_true(res);
    res = testFixture.OptiClusterGetsTheCorrectTag("opti_mcc12");
    expect_false(res);
  }
  test_that("Opticluster gets the correct number of bins") {
    OptiClusterTestFixture testFixture;
    bool res = testFixture.OptiClusterGetsTheCorrectNumberOfBins(1);
    expect_true(res);
    res = testFixture.OptiClusterGetsTheCorrectNumberOfBins(0);
    expect_false(res);
  }

}

#endif