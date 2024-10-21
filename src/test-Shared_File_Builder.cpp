//
// Created by Gregory Johnson on 8/30/24.
//
#include "TestHelpers/TestHelper.h"
#if DEBUG_RCPP
#include <testthat.h>
#include "Tests/SharedFileBuilderTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("SharedFile Builder") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("SharedFile Builder returns a SharedFile Object Successfully") {
        SharedFileBuilderTestFixture test;
        CountTableAdapter countTable;
        ListVector vector;
        vector.push_back("1");
        vector.push_back("2");
        vector.push_back("3");
        vector.push_back("4");
        vector.push_back("5");
        const std::vector<std::string> compounds{"1", "2", "3", "4", "5"};
        const std::vector<double> total{10, 20, 30, 40, 50};
        Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
            Rcpp::Named("Representative Sequence") = compounds,
            Rcpp::Named("total") = total,
            Rcpp::Named("nogroup") = total);
        countTable.CreateDataFrameMap(dataframe);
        bool res = test.TestBuildSharedFile(vector,countTable, true);
        expect_true(res);
        res = test.TestBuildSharedFile(vector,countTable, false);
        expect_false(res);
    }
}

#endif

