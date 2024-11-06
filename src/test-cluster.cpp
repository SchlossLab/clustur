//
// Created by gregj on 7/9/2024.
//

#include <testthat.h>
#include "Adapters/CountTableAdapter.h"
#include "Tests/ClusterTestFixture.h"
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/WeightedLinkage.h"
#include "MothurDependencies/AverageLinkage.h"
#include "MothurDependencies/SingleLinkage.h"
#include "MothurDependencies/CompleteLinkage.h"

context("Cluster algorithms") {
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("ClusterBins succesfully clusters bins") {
    ClusterTestFixture test_fixture;
      const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
      const std::vector<double> total{10, 20, 30, 40, 50, 60};
      const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
          Rcpp::Named("Representative Sequence") = compounds,
          Rcpp::Named("total") = total,
          Rcpp::Named("nogroup") = total);
    CountTableAdapter countTable;
    countTable.CreateDataFrameMap(dataframe);
    MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false, countTable);
    auto dMatrix = adapter.CreateSparseMatrix();
    ListVector listVector = adapter.CreateListVector();
    bool result = test_fixture.TestClusterNames(&listVector, &dMatrix, true);
    expect_true(result);
    result = test_fixture.TestClusterNames(&listVector, &dMatrix, false);
    expect_false(result);

  }
  test_that("Update map function does and does not fail") {
     ClusterTestFixture test_fixture;
      const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
      const std::vector<double> total{10, 20, 30, 40, 50, 60};
      const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
          Rcpp::Named("Representative Sequence") = compounds,
          Rcpp::Named("total") = total,
          Rcpp::Named("nogroup") = total);
      CountTableAdapter countTable;
      countTable.CreateDataFrameMap(dataframe);
    MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false, countTable);
    auto dMatrix = adapter.CreateSparseMatrix();
    ListVector listVector = adapter.CreateListVector();
    bool result = test_fixture.TestUpdateMap(&listVector, &dMatrix, true);
    expect_true(result);
    result = test_fixture.TestUpdateMap(&listVector, &dMatrix, false);
    expect_false(result);
  }
  test_that("Get Tag returns the tag that corresponds to the clustering type") {
    ClusterTestFixture test_fixture;
      const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
      const std::vector<double> total{10, 20, 30, 40, 50, 60};
      const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
          Rcpp::Named("Representative Sequence") = compounds,
          Rcpp::Named("total") = total,
          Rcpp::Named("nogroup") = total);
      CountTableAdapter countTable;
      countTable.CreateDataFrameMap(dataframe);
    MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false, countTable);
    auto dMatrix = adapter.CreateSparseMatrix();
    ListVector listVector = adapter.CreateListVector();

    auto* clustComplete = new CompleteLinkage(new RAbundVector(), &listVector,
        &dMatrix, 0.2, "furthest", -1);
    bool result = test_fixture.TestGetTagGetter(clustComplete, "fn");
    expect_true(result);
    result = test_fixture.TestGetTagGetter(clustComplete, "nn");
    expect_false(result);


    auto* clustSingle = new SingleLinkage(new RAbundVector(), &listVector,
        &dMatrix, 0.2, "nearest", -1);
    result = test_fixture.TestGetTagGetter(clustSingle, "nn");
    expect_true(result);
    result = test_fixture.TestGetTagGetter(clustSingle, "an");
    expect_false(result);

    auto* clustAverage = new AverageLinkage(new RAbundVector(), &listVector,
        &dMatrix, 0.2, "average", -1);
    result = test_fixture.TestGetTagGetter(clustAverage, "an");
    expect_true(result);
    result = test_fixture.TestGetTagGetter(clustAverage, "wn");
    expect_false(result);


    auto* clustWeighted = new WeightedLinkage(new RAbundVector(), &listVector,
      &dMatrix, 0.2, "weighted", -1);
    result = test_fixture.TestGetTagGetter(clustWeighted, "wn");
    expect_true(result);
    result = test_fixture.TestGetTagGetter(clustWeighted, "nn");
    expect_false(result);
    delete(clustWeighted);
    delete(clustComplete);
    delete(clustAverage);
    delete(clustSingle);
    }
    test_that("Get Tag returns the tag that corresponds to the clustering type") {
        ClusterTestFixture test_fixture;
      const std::vector<std::string> compounds{"1", "2", "3", "4", "5", "6"};
      const std::vector<double> total{10, 20, 30, 40, 50, 60};
      const Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
          Rcpp::Named("Representative Sequence") = compounds,
          Rcpp::Named("total") = total,
          Rcpp::Named("nogroup") = total);
      CountTableAdapter countTable;
      countTable.CreateDataFrameMap(dataframe);
      MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2, false, countTable);
      auto dMatrix = adapter.CreateSparseMatrix();
      ListVector listVector = adapter.CreateListVector();
      auto* clust = new SingleLinkage(new RAbundVector(), &listVector,
          &dMatrix, 0.2, "nearest", -1);

      bool result = test_fixture.TestGetSeqToBin(clust, std::map<std::string, int>());
      expect_true(result);
      result = test_fixture.TestGetSeqToBin(clust, std::map<std::string, int>{{"", 1}});
      expect_false(result);
      delete(clust);

    }
}
