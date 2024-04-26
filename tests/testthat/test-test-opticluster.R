test_that("Clustering returns proper results", {
  expected_df <- readRDS(here::here("tests/extdata/df_test_file.RDS"))
  matrix <- readRDS(here::here("tests/extdata/matrix_data.RDS"))
  df <- Opticluster::OptiCluster(matrix, 0.2, 2)
  expect_equal(class(df), "data.frame")
  expect_equal(df, expected_df)
})
