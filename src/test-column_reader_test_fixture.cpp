//
// Created by Gregory Johnson on 9/16/24.
//

#include <testthat.h>
#include "Tests/ColumnReaderTestFixture.h"
#include "RowData.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test contexta. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("ColumnDistanceMatrixReader Test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Column Reader reads column files") {
        ColumnReaderTestFixture fixture;
        bool result = fixture.TestReadColumnFile("", {}, true);
        expect_false(result);
    }
    test_that("Column Reader reads phylip files to rowdata information") {
        ColumnReaderTestFixture fixture;
        bool result = fixture.TestReadColumnFileToRowData("", {}, {});
        expect_false(result);
    }

    test_that("Column Reader properly creates a sparse distance matrix"){
        ColumnReaderTestFixture fixture;
        bool result = fixture.TestGetDistanceMatrix("", {}, false);
        expect_true(result);

    }
      test_that("Column Reader properly creates a list vector"){
        ColumnReaderTestFixture fixture;
        bool result = fixture.TestGetListVector("", {}, false);
        expect_true(result);

    }
}
