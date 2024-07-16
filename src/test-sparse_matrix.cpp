//
// Created by gregj on 7/9/2024.
//


#include <testthat.h>
#include "Tests/SparseMatrixTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("Spase Distance Matrix") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("The clear function removes data in the sequence vector") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestClear(true);
        expect_true(result);
        result = fixture.TestClear(false);
        expect_false(result);
    }
    test_that("The resize function, resize the data") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestResize(4, 4);
        expect_true(result);
        result = fixture.TestResize(4, 3);
        expect_false(result);
        result = fixture.TestResize(0, 0);
        expect_false(result);
    }
    test_that("Prints returns whether it printed the data or not") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestPrint(false, true);
        expect_true(result);
        result = fixture.TestPrint(true, false);
        expect_true(result);
        result = fixture.TestPrint(true, true);
        expect_false(result);
    }
     test_that("We can properly add cells to the sparse distance matrix") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestAddCell(2, PDistCell(1,2.0), true);
        expect_true(result);
        result = fixture.TestAddCell(2, PDistCell(1,2.0), false);
        expect_false(result);
    }
     test_that("Heap Comparator compares two PDistCells accurately") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestHeapComparator(PDistCell(2,2), PDistCell(1,3), true);
        expect_true(result);
        result = fixture.TestHeapComparator(PDistCell(2,4), PDistCell(1,3), true);
        expect_false(result);
    }
     test_that("RmCells remove cells correctly") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestRemoveCell(1,1,true);
        expect_true(result);
    }
    test_that("Test Adding Cells sorted works") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestAddCellSorted(2, PDistCell(1,2.0), true);
        expect_true(result);
        result = fixture.TestAddCellSorted(2, PDistCell(1,2.0), false);
        expect_false(result);
    }
    test_that("Test Getting the smallest cells gets the smallest cell value") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestGetSmallestCell(1, 0);
        expect_true(result);
        result = fixture.TestGetSmallestCell(1, 1);
        expect_false(result);
    }
    test_that("Test Getting the smallest distance returns the smallest distance") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestGetSmallestDist(0);
        expect_true(result);
        result = fixture.TestGetSmallestDist(1);
        expect_false(result);
    }
    test_that("Test Updating the Cell Compliment") {
        SparseMatrixTestFixture fixture;
        bool result = fixture.TestUpdateCellCompliment(2,2, true);
        expect_true(result);
        result = fixture.TestUpdateCellCompliment(2,2, false);
        expect_false(result);
    }

}