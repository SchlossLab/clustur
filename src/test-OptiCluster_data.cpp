//
// Created by Gregory Johnson on 8/30/24.
//
#include "TestHelpers/TestHelper.h"
#if DEBUG_RCPP
#include <testthat.h>
#include "Tests/OptiClusterDataTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("OptiCluster Data") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Cluster Data Adds to list properly") {
        OptiClusterDataTestFixture testFixture;
        OptiClusterInformation information{3, "as", "0.011"};
        bool res = testFixture.CheckAddToData(information, true);
        expect_true(res);
        OptiClusterInformation information2{3, "as", ""};
        testFixture.CheckAddToData(information2, false);
        expect_true(res);
    }
}



#endif
