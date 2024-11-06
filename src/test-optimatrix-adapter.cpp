

#include <testthat.h>
#include "Tests/OptimatrixAdapterTestFixture.h"
#include "MothurDependencies/ColumnDistanceMatrixReader.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
//  OptiMatrix* ConvertToOptimatrix(const std::vector<int>
//        &xPosition, const std::vector<int>& yPosition, const std::vector<double>& data, int rowSize, int colSize);
// void Check()
// {
//   OptimatrixAdapterTestFixture fixture;
//   std::vector<int> xVals = std::vector<int>{1,1,1,2,2,3,5};
//   std::vector<int> yVals = std::vector<int>{2,3,5,3,5,5,5};
//   std::vector<double> data = std::vector<double>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
//   bool result = fixture.TestOptimatrixClosenessReturnsCorrectValue(xVals, yVals, data, 4);
// }
context("Optimatrix Adapter Test") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("OptimatrixAdapter returns a optimatrix") {
    OptimatrixAdapterTestFixture fixture;
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
    Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
    const Rcpp::Function test_path = pkg["test_path"];
    const Rcpp::Function read_count = clustur["read_count"];
    const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_column.dist"));
    const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
    const Rcpp::DataFrame df = read_count(countTablePath);
    ColumnDistanceMatrixReader reader(0.2, false);
    reader.CreateCountTableAdapter(df);
    reader.Read(path);
    bool result = fixture.TestOptimatrixReturnsNotNullValues(reader.GetSparseMatrix(), reader.GetListVector());
    expect_true(result);
  }
  test_that("OptimatrixAdapter returns the correct number of Closeness Values"){
    OptimatrixAdapterTestFixture fixture;
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
    Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
    const Rcpp::Function test_path = pkg["test_path"];
    const Rcpp::Function read_count = clustur["read_count"];
    const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_column.dist"));
    const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
    const Rcpp::DataFrame df = read_count(countTablePath);
    ColumnDistanceMatrixReader reader(0.2, false);
    reader.CreateCountTableAdapter(df);
    reader.Read(path);
    bool result = fixture.TestOptimatrixClosenessReturnsCorrectValue(reader.GetSparseMatrix(), reader.GetListVector(), 86);
    expect_true(result);
  }
  test_that("OptimatrixAdapter returns the correct number of Singletons"){
     OptimatrixAdapterTestFixture fixture;
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
    Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
    const Rcpp::Function test_path = pkg["test_path"];
    const Rcpp::Function read_count = clustur["read_count"];
    const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_column.dist"));
    const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
    const Rcpp::DataFrame df = read_count(countTablePath);
    ColumnDistanceMatrixReader reader(0.2, false);
    reader.CreateCountTableAdapter(df);
    reader.Read(path);
    bool result = fixture.TestOptimatrixSingletonReturnsCorrectValue(reader.GetSparseMatrix(), reader.GetListVector(), 12);
    expect_true(result);
  }
  test_that("OptimatrixAdapter returns the correct number of Names"){
    OptimatrixAdapterTestFixture fixture;
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
    Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
    const Rcpp::Function test_path = pkg["test_path"];
    const Rcpp::Function read_count = clustur["read_count"];
    const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_column.dist"));
    const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
    const Rcpp::DataFrame df = read_count(countTablePath);
    ColumnDistanceMatrixReader reader(0.2, false);
    reader.CreateCountTableAdapter(df);
    reader.Read(path);
    bool result = fixture.TestOptimatrixNameListReturnsCorrectValue(reader.GetSparseMatrix(), reader.GetListVector(), 98);
    expect_true(result);
  }

}

