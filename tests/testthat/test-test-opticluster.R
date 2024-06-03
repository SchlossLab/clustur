test_that("Clustering returns proper results", {
  expected_df <- readRDS(test_path("extdata","df_test_file.RDS"))
  matrix <- readRDS(test_path("extdata","matrix_data.RDS"))
  df <- Opticluster::opti_cluster(matrix, 0.2, 2)
  expect_equal(class(df), "data.frame")
  # expect_equal(df, expected_df)
})