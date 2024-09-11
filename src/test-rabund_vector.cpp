//
// Created by Gregory Johnson on 7/12/24.
//
#include <testthat.h>
#include "Tests/RAbundVectorTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("RAdbundVector Test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("RAbund vector can properly clear data") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestClear(true);
        expect_true(result);
        result = fixture.TestClear(false);
        expect_false(result);
    }
        test_that("RAbund vector can get abundance given index") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestGet(1, 1);
        expect_true(result);
        result = fixture.TestGet(1, 0);
        expect_false(result);
    }
        test_that("RAbund vector can quicksort") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestQuicksort(true);
        expect_true(result);
        result = fixture.TestQuicksort(false);
        expect_false(result);
    }
        test_that("RAbund vector can properly remove data given a bin") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestRemove(1,true);
        expect_true(result);
        result = fixture.TestRemove(1,false);
        expect_false(result);
    }
        test_that("RAbund vector can properly resize its data") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestResize(3,true);
        expect_true(result);
        result = fixture.TestResize(3,false);
        expect_false(result);
    }
        test_that("RAbund vector can set new bin sizes inside of its set function") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestSet(1,2, true);
        expect_true(result);
        result = fixture.TestSet(1,2, false);
        expect_false(result);
    }    test_that("RAbund vector can properly return its correct size") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestSize(5);
        expect_true(result);
        result = fixture.TestSize(4);
        expect_false(result);
    }
        test_that("RAbund vector can properly return its data") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestGetOverload({1,1,1,1,1});
        expect_true(result);
        result = fixture.TestGetOverload({1,1,1,1});
        expect_false(result);
    }
       test_that("RAbund vector push_back will add to the vector properly") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestPushBack(2, true);
        expect_true(result);
        result = fixture.TestPushBack(2, false);
        expect_false(result);
    }
       test_that("RAbund vector can properly return its max rank") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestGetMaxRank(1);
        expect_true(result);
        result = fixture.TestGetMaxRank(3);
        expect_false(result);
    }
       test_that("RAbund vector can properly return its number of bins") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestGetNumBins(5);
        expect_true(result);
        result = fixture.TestGetNumBins(4);
        expect_false(result);
    }
           test_that("RAbund vector can properly return its number of sequences") {
        RAbundVectorTestFixture fixture;
        bool result = fixture.TestGetNumSeqs(5);
        expect_true(result);
        result = fixture.TestGetNumSeqs(1);
        expect_false(result);
    }

}