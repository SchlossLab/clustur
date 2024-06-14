//
// Created by Gregory Johnson on 5/7/24.
//
//
// Created by Gregory Johnson on 5/2/24.
//
//
// Created by Gregory Johnson on 5/1/24.
//

#include "TestHelpers/TestHelper.h"
#include "Tests/OptiDataTestFixture.h"
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

context("Optidata test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Optidata gets the correct close sequences") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestGetCloseSequencesReturnCorrectData(1, {0,2,4});
        expect_true(result);
        result = testFixture.TestGetCloseSequencesReturnCorrectData(1, {0,1,4});
        expect_false(result);
    }
    test_that("Optidata Correct test for close data") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestIsCloseReturnsIfDataIsClose(0, 1, true);
        expect_true(result);
        result = testFixture.TestIsCloseReturnsIfDataIsClose(0, 2, true);
        expect_false(result);
    }
    test_that("Test GetNumsClose Returns the correct number of close variables") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestGetNumCloseReturnsCorrectData(1, 3);
        expect_true(result);
        result = testFixture.TestGetNumCloseReturnsCorrectData(1, 2);
        expect_false(result);

    }
    test_that("Optidata Returns the correct sequence name") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestGetNameReturnsTheProperName(1, "1");
        expect_true(result);
        result = testFixture.TestGetNameReturnsTheProperName(1, "3");
        expect_false(result);

    }
    test_that("Optidata Returns the correct number of sequences") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestGetNumSeqsReturnsCorrectData(4);
        expect_true(result);
        result = testFixture.TestGetNumSeqsReturnsCorrectData(3);
        expect_false(result);
    }
    test_that("Optidata Returns the correct number of singletons") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestGetNumSingletonsReturnsCorrectData(1);
        expect_true(result);
        result = testFixture.TestGetNumSingletonsReturnsCorrectData(2);
        expect_false(result);
    }
    test_that("Optidata GetList returns a proper listVector and the correct number of numSeqs") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestGetListSingleReturnsCorrectData(1);
        expect_true(result);
        result = testFixture.TestGetListSingleReturnsCorrectData(2);
        expect_false(result);
    }
    test_that("Optidata Returns Data for ref sequences") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestGetCloseRefSeqsReturnsCorrectData(1, std::set<long long>());
        expect_true(result);
        result = testFixture.TestGetCloseRefSeqsReturnsCorrectData(1, std::set<long long>{1,2});
        expect_false(result);
    }
    test_that("Optidata IsCloseFit returns false") {
        OptiDataTestFixture testFixture;
        bool result = testFixture.TestIsCloseFitReturnsCorrectData(1, 1, false, false);
        expect_true(result);
        result = testFixture.TestIsCloseFitReturnsCorrectData(1, 1, false, true);
        expect_false(result);
    }



}

#endif