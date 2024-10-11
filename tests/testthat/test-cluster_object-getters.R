test_that("get_cutoff returns the proper cutoff", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, method = "furthest")
  expect_true(get_cutoff(df_clusters) == 0.19)
  expect_false(get_cutoff(df_clusters) == 0.4)
})

test_that("get_clusters returns the cluster dataframe", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, method = "furthest")
  expect_true(class(get_clusters(df_clusters)) == "data.frame")
  expect_true(all(names(get_clusters(df_clusters))
                  == c("label", "bins", "otu")))
})

test_that("get_metrics returns the metrics for opticlustered clusters", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, method = "furthest")
  expect_error(get_metrics(df_clusters))
  df_clusters <- cluster(distance_data, method = "opti")
  metrics <- get_metrics(df_clusters)
  expect_true(length(metrics) == 2)
  expect_true(ncol(metrics$metrics) == 14)
  expect_true(ncol(metrics$other_metrics) == 17)
})

test_that("get_shared returns the shared dataframe", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df_clusters <- cluster(distance_data, method = "furthest")
  expect_true(names(count_table)[3:length(count_table)] %in%
                get_shared(df_clusters)$samples)
})
