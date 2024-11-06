//
// Created by gregj on 7/9/2024.
//

#include <testthat.h>
#include "Tests/PhylipReaderTestFixture.h"
#include "Adapters/MatrixAdapter.h"
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
        Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
        const Rcpp::Function test_path = pkg["test_path"];
        const Rcpp::Function read_count = clustur["read_count"];
        const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_phylip.dist"));
        const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
        const Rcpp::DataFrame df = read_count(countTablePath);
        ReadPhylipMatrix reader(0.2, false);
        bool result = fixture.TestReadPhylipFile(df, path, true);
        expect_true(result);
        result = fixture.TestReadPhylipFile(df, "", true);
        expect_false(result);
    }
    test_that("Phylip Reader properly creates a sparse distance matrix"){
        PhylipReaderTestFixture fixture;
        Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
        Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
        const Rcpp::Function test_path = pkg["test_path"];
        const Rcpp::Function read_count = clustur["read_count"];
        const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_phylip.dist"));
        const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
        const Rcpp::DataFrame df = read_count(countTablePath);
        ReadPhylipMatrix reader(0.2, false);
        bool result = fixture.TestGetSparseMatrix(df, path, true);
        expect_true(result);
        result = fixture.TestGetSparseMatrix(df, path, false);
        expect_false(result);

    }
      test_that("Phylip Reader properly creates a list vector"){
        PhylipReaderTestFixture fixture;
        Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
        Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
        const Rcpp::Function test_path = pkg["test_path"];
        const Rcpp::Function read_count = clustur["read_count"];
        const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_phylip.dist"));
        const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
        const Rcpp::DataFrame df = read_count(countTablePath);
        ReadPhylipMatrix reader(0.2, false);
        bool result = fixture.TestGetListVector(df, path, 98);
        expect_true(result);

        result = fixture.TestGetListVector(df, path, 3);
        expect_false(result);

    }
}
