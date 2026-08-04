test_that("get_label returns the proper cutoff", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, cutoff, method = "furthest")
  expect_true(get_cutoff(df_clusters) == 0.19)
  expect_false(get_cutoff(df_clusters) == 0.4)

  # Error Checking
  expect_error(get_cutoff(""))
  expect_error(get_cutoff(list()))
})

test_that("get_clusters returns the cluster dataframe", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, cutoff, method = "furthest")
  expect_true(class(get_bins(df_clusters)) == "data.frame")
  expect_true(all(names(get_bins(df_clusters))
                  %in%  c("feature", "bin")))

  # Error Checking
  expect_error(get_bins(""))
  expect_error(get_bins(list()))
})

test_that("get_metrics returns the metrics for opticlustered clusters", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, cutoff, method = "furthest")
  expect_error(get_metrics(df_clusters))
  df_clusters <- cluster(distance_data, cutoff, method = "opticlust")
  metrics <- get_metrics(df_clusters)
  expect_true(length(metrics) == 2)
  expect_true(ncol(metrics$metrics) == 14)
  expect_true(ncol(metrics$iteration_metrics) == 17)

  # Error Checking
  expect_error(get_metrics(""))
  expect_error(get_metrics(list()))
})

test_that("get_shared returns the shared dataframe", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, cutoff, method = "furthest")
  expect_true(names(count_table)[3:length(count_table)] %in%
                get_abundance(df_clusters)$samples)

  # Error Checking
  expect_error(get_shared(""))
  expect_error(get_shared(list()))
})
