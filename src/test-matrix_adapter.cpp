//
// Created by Gregory Johnson on 7/12/24.
//

#include <testthat.h>
#include "Tests/MatrixAdapterTestFixture.h"
#include "Adapters/MatrixAdapter.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("MatrixAdapter Test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.

     test_that("Matrix Adapter can get its list vector") {
         const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
         const std::vector<double> total{10, 20, 30, 40, 50, 60};
         const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
             Rcpp::Named("Representative Sequence") = compounds,
             Rcpp::Named("total") = total,
             Rcpp::Named("nogroup") = total);
         CountTableAdapter countTable;
         countTable.CreateDataFrameMap(dataframe);

         const MatrixAdapter adapter = MatrixAdapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25},
             0.2, false, countTable);

         MatrixAdapterTestFixture fixture;
         bool result = fixture.TestGetListVector(adapter, true);
         expect_true(result);
         result = fixture.TestGetListVector(adapter, false);
         expect_false(result);
     }
     test_that("Matrix Adapter can create sparse matrices") {
         const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
         const std::vector<double> total{10, 20, 30, 40, 50, 60};
         const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
             Rcpp::Named("Representative Sequence") = compounds,
             Rcpp::Named("total") = total,
             Rcpp::Named("nogroup") = total);
         CountTableAdapter countTable;
         countTable.CreateDataFrameMap(dataframe);

         const MatrixAdapter adapter = MatrixAdapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25},
             0.2, false, countTable);

         MatrixAdapterTestFixture fixture;
         bool result = fixture.TestCreateSparseMatrix(adapter, 6);
         expect_true(result);
         result = fixture.TestCreateSparseMatrix(adapter, 0);
         expect_false(result);
     }
 }
