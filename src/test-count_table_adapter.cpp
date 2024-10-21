//
// Created by Gregory Johnson on 8/30/24.
//
#include "TestHelpers/TestHelper.h"
#if DEBUG_RCPP
#include <testthat.h>
#include "Tests/CountTableAdapterTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("CountTable") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Count table properly creates a dataframe map") {
        CountTableAdapterTestFixture testFixture;
        const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
        const std::vector<double> total{10, 20, 30, 40, 50, 60};
        const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
        Rcpp::Named("Representative Sequence") = compounds,
        Rcpp::Named("total") = total,
        Rcpp::Named("nogroup") = total);
        bool res = testFixture.TestCreateDataFrameMap(dataframe, true);
        expect_true(res);
        const Rcpp::DataFrame dataframe2 = Rcpp::DataFrame::create(
        Rcpp::Named("Representative Sequence") = compounds,
        Rcpp::Named("total") = total);
        res = testFixture.TestCreateDataFrameMap(dataframe2, false);
        expect_true(res);
    }
    test_that("Count table properly finds abundances") {
        CountTableAdapterTestFixture testFixture;
        bool res = testFixture.TestFindAbundance("2", true);
        expect_true(res);
        res = testFixture.TestFindAbundance("10", false);
        expect_true(res);
    }
    test_that("Count table properly finds abundance based on group") {
        CountTableAdapterTestFixture testFixture;
        bool res = testFixture.TestFindAbundanceBasedOnGroup("nogroup", "2", 20);
        expect_true(res);
        res = testFixture.TestFindAbundanceBasedOnGroup("nogroup", "2", 30);
        expect_false(res);
    }
    test_that("Count table properly returns column by name") {
        CountTableAdapterTestFixture testFixture;
        bool res = testFixture.TestGetColumnByName("nogroup", true);
        expect_true(res);
        res = testFixture.TestGetColumnByName("nogroup1", false);
        expect_true(res);
    }
    test_that("Count table properly returns groups") {
        CountTableAdapterTestFixture testFixture;
        bool res = testFixture.TestGetGroups({"nogroup"});
        expect_true(res);
        res = testFixture.TestGetGroups({});
        expect_false(res);
        res = testFixture.TestGetGroups({"1"});
        expect_false(res);
    }
    test_that("Count table properly returns the correct name given the index") {
        CountTableAdapterTestFixture testFixture;
        bool res = testFixture.TestGetNameByIndex(1, "2");
        expect_true(res);
        res = testFixture.TestGetNameByIndex(7, "");
        expect_true(res);
        res = testFixture.TestGetNameByIndex(2, "");
        expect_false(res);
    }


}



#endif
