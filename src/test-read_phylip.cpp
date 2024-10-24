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

// Initialize a unit test contexta. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("PhylipReader Test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Phylip Reader reads phylip files") {
        PhylipReaderTestFixture fixture;
        Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
        const Rcpp::Function func = pkg["test_path"];
        const std::string path = Rcpp::as<std::string>(func("extdata", "amazon_phylip.dist"));
        ReadPhylipMatrix reader(0.2, false);
        reader.Read(path);
        bool result = fixture.TestReadPhylipFile(path, reader.GetListVector()->getNumSeqs());
        expect_true(result);
        result = fixture.TestReadPhylipFile("", 1);
        expect_false(result);
    }
    test_that("Phylip Reader reads phylip files to rowdata information") {
        PhylipReaderTestFixture fixture;
        Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
        const Rcpp::Function func = pkg["test_path"];
        const std::string path = Rcpp::as<std::string>(func("extdata", "amazon_phylip.dist"));
        ReadPhylipMatrix reader(0.2, false);
        bool result = fixture.TestReadPhylipFileToRowData(path, reader.ReadToRowData(path));
        expect_true(result);
        result = fixture.TestReadPhylipFileToRowData("", std::vector<RowData>(1));
        expect_false(result);
    }
    test_that("Phylip Reader reads phylip file from created RowData Vectors") {
        PhylipReaderTestFixture fixture;
        const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
        const std::vector<double> total{10, 20, 30, 40, 50, 60};
        const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
            Rcpp::Named("Representative Sequence") = compounds,
            Rcpp::Named("total") = total,
            Rcpp::Named("nogroup") = total);
        CountTableAdapter countTable;
        countTable.CreateDataFrameMap(dataframe);
        MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false, countTable);
        std::vector<RowData> squareMatrix = adapter.DistanceMatrixToSquareMatrix();
        bool result = fixture.TestReadPhylipFileFromRowData(squareMatrix, 6);
        expect_true(result);
        result = fixture.TestReadPhylipFileFromRowData(squareMatrix, 3);
        expect_false(result);
    }
    test_that("Phylip Reader properly creates a sparse distance matrix"){
        PhylipReaderTestFixture fixture;
        const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
        const std::vector<double> total{10, 20, 30, 40, 50, 60};
        const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
            Rcpp::Named("Representative Sequence") = compounds,
            Rcpp::Named("total") = total,
            Rcpp::Named("nogroup") = total);
        CountTableAdapter countTable;
        countTable.CreateDataFrameMap(dataframe);
        MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false, countTable);
        std::vector<RowData> squareMatrix = adapter.DistanceMatrixToSquareMatrix();
        bool result = fixture.TestGetDistanceMatrix(squareMatrix, true);
        expect_true(result);

        MatrixAdapter adapterTwo({}, {}, {}, 0.2, false, countTable);
        std::vector<RowData> squareMatrixTwo = adapterTwo.DistanceMatrixToSquareMatrix();
        result = fixture.TestGetDistanceMatrix(squareMatrixTwo, true);
        expect_false(result);

    }
      test_that("Phylip Reader properly creates a list vector"){
        PhylipReaderTestFixture fixture;
        const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
        const std::vector<double> total{10, 20, 30, 40, 50, 60};
        const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
            Rcpp::Named("Representative Sequence") = compounds,
            Rcpp::Named("total") = total,
            Rcpp::Named("nogroup") = total);
        CountTableAdapter countTable;
        countTable.CreateDataFrameMap(dataframe);
        MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false, countTable);
        std::vector<RowData> squareMatrix = adapter.DistanceMatrixToSquareMatrix();
        bool result = fixture.TestGetListVector(adapter.DistanceMatrixToSquareMatrix(), 6);
        expect_true(result);

        MatrixAdapter adapterTwo({}, {}, {}, 0.2, false, countTable);
        std::vector<RowData> squareMatrixTwo = adapterTwo.DistanceMatrixToSquareMatrix();
        result = fixture.TestGetListVector(squareMatrixTwo, 3);
        expect_false(result);

    }
}
