test_that("Clustering returns proper results", {
  expected_df <- readRDS(test_path("extdata","df_test_file.RDS"))
  matrix <- readRDS(test_path("extdata","matrix_data.RDS"))
  df <- Opticluster::opti_cluster(matrix, 0.2, 2, FALSE)
  df$cluster$exists <- do.call(paste0, df$cluster) %in% do.call(paste0, expected_df)
  expect_equal(class(df$cluster), "data.frame")
  expect_equal(class(df$cluster_metrics), "data.frame")
  expect_equal(class(df$other_cluster_metrics), "data.frame")
  expect_true(all(df$exists == TRUE))
})

test_that("Normal Cluster is able to properly cluster data",
{
  expected_df <- readRDS(test_path("extdata","df_test_file.RDS"))
  matrix <- readRDS(test_path("extdata","matrix_data.RDS"))
  
  cluster_furthest <- cluster(matrix, 0.2, "furthest")
  cluster_average <- cluster(matrix, 0.2, "average")
  cluster_weighted <- cluster(matrix, 0.2, "weighted")
  cluster_nearest <- cluster(matrix, 0.2, "nearest")

  expect_true(any(class(cluster_furthest) == "character"))
  expect_true(any(class(cluster_average) == "character"))
  expect_true(any(class(cluster_weighted) == "character"))
  expect_true(any(class(cluster_nearest) == "character"))

})