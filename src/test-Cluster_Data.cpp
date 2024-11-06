//
// Created by Gregory Johnson on 8/30/24.
//
#include <testthat.h>
#include "Tests/ClusterDataTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("Cluster Data") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Cluster Data Adds to list properly") {
        ClusterDataTestFixture test;
        ClusterInformation info{"0.011", 3, "qa"};
        bool res = test.CheckAddToData(info, true);
        expect_true(res);
        ClusterInformation info2 {"", 3, "qa"};
        res = test.CheckAddToData(info2, false);
        expect_true(res);
    }
}
