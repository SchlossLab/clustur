//
// Created by Gregory Johnson on 8/30/24.
//
#include "TestHelpers/TestHelper.h"
#if DEBUG_RCPP
#include <testthat.h>
#include "Tests/ClusterExportTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("Cluster Export") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Cluster Exports sets list vector properly") {
        ClusterExportTestFixture testFixture;
        const auto vec = new ListVector();
        const std::string label = "0.01";
        bool res = testFixture.TestSetListVector(*vec, label, true);
        expect_true(res);
        const auto vec2 = new ListVector();
        res = testFixture.TestSetListVector(*vec2, label, false);
        expect_false(res);
    }
    test_that("Cluster Exports returns list vector properly") {
        ClusterExportTestFixture testFixture;
        bool res = testFixture.TestGetListVector(true);
        expect_true(res);
        res = testFixture.TestGetListVector(false);
        expect_false(res);
    }


}



#endif