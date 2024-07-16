//
// Created by gregj on 7/9/2024.
//

#include <testthat.h>
#include "Tests/PhylipReaderTestFixture.h"
#include "Adapters/MatrixAdapter.h"
#include "RowData.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("PhylipReader Test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Phylip Reader reads phylip files") {
        PhylipReaderTestFixture fixture;
        //bool result = fixture.TestReadPhylipFile
        //("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/testthat/extdata/matrix_phylip.txt", true);
        // expect_true(result);
        bool result = fixture.TestReadPhylipFile("", true);
        expect_false(result);
    }
    test_that("Phylip Reader reads phylip file from created RowData Vectors") {
        PhylipReaderTestFixture fixture;
        MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2);
        std::vector<RowData> squareMatrix = adapter.DistanceMatrixToSquareMatrix();
        bool result = fixture.TestReadPhylipFile(squareMatrix, true);
        expect_true(result);
        result = fixture.TestReadPhylipFile(squareMatrix, false);
        expect_false(result);
    }
    test_that("Phylip Reader properly creates a sparse distance matrix"){
        PhylipReaderTestFixture fixture;
        MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2);
        std::vector<RowData> squareMatrix = adapter.DistanceMatrixToSquareMatrix();
        bool result = fixture.TestGetDistanceMatrix(squareMatrix, true);
        expect_true(result);

        MatrixAdapter adapterTwo({}, {}, {}, 0.2);
        std::vector<RowData> squareMatrixTwo = adapterTwo.DistanceMatrixToSquareMatrix();
        result = fixture.TestGetDistanceMatrix(squareMatrixTwo, true);
        expect_false(result);

    }
      test_that("Phylip Reader properly creates a list vector"){
        PhylipReaderTestFixture fixture;
        MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2);
        std::vector<RowData> squareMatrix = adapter.DistanceMatrixToSquareMatrix();
        bool result = fixture.TestGetListVector(adapter.DistanceMatrixToSquareMatrix(), true);
        expect_true(result);

        MatrixAdapter adapterTwo({}, {}, {}, 0.2);
        std::vector<RowData> squareMatrixTwo = adapterTwo.DistanceMatrixToSquareMatrix();
        result = fixture.TestGetListVector(squareMatrixTwo, true);
        expect_false(result);

    }
}