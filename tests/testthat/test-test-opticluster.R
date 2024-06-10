test_that("Clustering returns proper results", {
  expected_df <- readRDS(test_path("extdata","df_test_file.RDS"))
  matrix <- readRDS(test_path("extdata","matrix_data.RDS"))
  df <- Opticluster::opti_cluster(matrix, 0.2, 2)
  df$exists <- do.call(paste0, df) %in% do.call(paste0, expected_df)
  expect_equal(class(df), "data.frame")
  expect_true(all(df$exists == TRUE))
})
