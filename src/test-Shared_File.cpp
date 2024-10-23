//
// Created by Gregory Johnson on 8/30/24.
//
#include "TestHelpers/TestHelper.h"
#if DEBUG_RCPP
#include <testthat.h>
#include "Tests/SharedFileTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("SharedFile") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("SharedFile builds out a proper dataframe") {
        SharedFileTestFixture fixture;
        SharedAbundance abundance{"grp1", "1", "0.000", 10};
        SharedAbundance abundance2{"grp2", "1", "0.000", 20};
        const std::vector<std::string> compounds{"1", "2", "3", "4", "5"};
        const std::vector<double> total{10, 20, 30, 40, 50};
        Rcpp::DataFrame dataFrame = Rcpp::DataFrame::create(Rcpp::Named("samples") = compounds,
                                                    Rcpp::Named("otu") = compounds,
                                                    Rcpp::Named("abundance") = total);

        bool res = fixture.TestSharedFilePrintData({abundance, abundance2}, dataFrame);
        expect_true(res);

        Rcpp::DataFrame dataFrame2 = Rcpp::DataFrame::create(Rcpp::Named("samples") = compounds,
                                                    Rcpp::Named("otu") = compounds,
                                                    Rcpp::Named("random") = total);
        res = fixture.TestSharedFilePrintData({abundance, abundance2}, dataFrame2);
        expect_false(res);
    }
}

#endif

