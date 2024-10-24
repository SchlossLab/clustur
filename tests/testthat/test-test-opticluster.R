test_that("opticluster returns four dataframes", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df <- cluster(distance_data, method = "opticlust")
  csv <- read.csv(test_path("extdata", "abundance_results_opticluster.csv"))
  expect_equal(class(df$cluster), "data.frame")
  expect_equal(class(df$cluster_metrics), "data.frame")
  expect_equal(class(df$iteration_metrics), "data.frame")
  expect_equal(class(df$abundance), "data.frame")
  expect_true(all(csv$label == df$label))
})

test_that("other clustering methods only return two dataframes", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)

  df_furthest <- cluster(distance_data, method = "furthest")
  df_nearest <- cluster(distance_data, method = "nearest")
  df_average <- cluster(distance_data, method = "average")
  df_weighted <- cluster(distance_data, method = "weighted")

  expect_equal(class(df_furthest$cluster), "data.frame")
  expect_equal(class(df_furthest$abundance), "data.frame")

  expect_equal(class(df_nearest$cluster), "data.frame")
  expect_equal(class(df_nearest$abundance), "data.frame")

  expect_equal(class(df_average$cluster), "data.frame")
  expect_equal(class(df_average$abundance), "data.frame")

  expect_equal(class(df_weighted$cluster), "data.frame")
  expect_equal(class(df_weighted$abundance), "data.frame")
})

test_that("Clustering fails when arguments are incorrect", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  expect_error(cluster(distance_data, "neighbor"))
  expect_error(cluster(cutoff, "opticlust"))
})

test_that("Opticluster works with phylip files", {

  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_phylip.dist"),
                             count_table, cutoff, FALSE)

  df_opti <- cluster(distance_data, method = "opticlust")

  expect_true(nrow(df_opti$cluster) == 29 && nrow(df_opti$abundance == 29))
  expect_true(length(df_opti) == 5)

})

test_that("Example Path returns the correct path", {
  expect_true(class(example_path("amazon_phylip.dist")) == "character")
  expect_error(example_path("a"))
  expect_true(length(example_path()) == 4)
})

test_that("cluster works via phylip file", {

  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_phylip.dist"),
                             count_table, cutoff, FALSE)

  df_furthest <- cluster(distance_data, method = "furthest")
  df_nearest <- cluster(distance_data, method = "nearest")
  df_average <- cluster(distance_data, method = "average")
  df_weighted <- cluster(distance_data, method = "weighted")

  expect_true(nrow(df_furthest$cluster) == 31
              && nrow(df_furthest$abundance == 31))
  expect_true(nrow(df_nearest$cluster) == 20
              && nrow(df_furthest$abundance == 20))
  expect_true(nrow(df_average$cluster) == 35
              && nrow(df_furthest$abundance == 35))
  expect_true(nrow(df_weighted$cluster) == 36
              && nrow(df_furthest$abundance == 36))
})

test_that("Opticluster works with column files", {

  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)

  df_opti <- cluster(distance_data, method = "opticlust")


  expect_true(nrow(df_opti$cluster) == 31 && nrow(df_opti$abundance == 31))
})

test_that("Normal cluster works via column file", {

  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)

  df_furthest <- cluster(distance_data, method = "furthest")
  df_nearest <- cluster(distance_data, method = "nearest")
  df_average <- cluster(distance_data, method = "average")
  df_weighted <- cluster(distance_data, method = "weighted")

  expect_true(nrow(df_furthest$cluster) == 33
              && nrow(df_furthest$abundance == 33))
  expect_true(nrow(df_nearest$cluster) == 22
              && nrow(df_furthest$abundance == 22))
  expect_true(nrow(df_average$cluster) == 37
              && nrow(df_furthest$abundance == 37))
  expect_true(nrow(df_weighted$cluster) == 38
              && nrow(df_furthest$abundance == 38))
})


test_that("Amazon Data from mothur clusters properly", {

  result <- readRDS(test_path("extdata", "amazon_opti_results.RDS"))
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, 0.2, FALSE)
  data <- cluster(distance_data, method = "opticlust")
  expect_true(nrow(data$cluster) == nrow(result$cluster))
})

test_that("Read count can read sparse and normal count tables", {
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  sparse_count_table <- read_count(test_path("extdata",
                                             "amazon_sparse.count_table"))

  expect_true("data.frame" %in% class(count_table))
  expect_true(nrow(count_table) == 98)
  expect_true(ncol(count_table) == 3)


  expect_true("data.frame" %in% class(sparse_count_table))
  expect_true(nrow(sparse_count_table) == 5)
  expect_true(ncol(sparse_count_table) == 12)
})

test_that("Read dist can read column, phylip files, and sparse matrices", {
  set.seed(123)
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data_column <- read_dist(test_path("extdata", "amazon_column.dist"),
                                    count_table, 0.2, FALSE)
  distance_data_phylip <- read_dist(test_path("extdata", "amazon_phylip.dist"),
                                    count_table, 0.2, FALSE)
  i_values <- as.integer(1:100)
  j_values <- as.integer(sample(1:100, 100, TRUE))
  x_values <- as.numeric(runif(100, 0, 1))
  s_matrix <- create_sparse_matrix(i_values, j_values, x_values)
  sparse_count <- data.frame(Representative_Sequence = 1:100,
                             total = rep(1, times = 100))
  distance_data_sparse <- read_dist(s_matrix, sparse_count, 0.2, FALSE)
  expect_error(read_dist("", count_table, 0.2, FALSE))
  expect_true(nrow(get_distance_df(distance_data_column)) == 962)
  expect_true(nrow(get_distance_df(distance_data_phylip)) == 1048)
  expect_true(nrow(get_distance_df(distance_data_sparse)) == 9827)
})


test_that("We can determine if a file is phylip or not", {
  is_not_phylip <-
    DetermineIfPhylipOrColumnFile(test_path("extdata", "amazon_column.dist"))
  is_phylip <-
    DetermineIfPhylipOrColumnFile(test_path("extdata", "amazon_phylip.dist"))
  expect_true(is_phylip)
  expect_false(is_not_phylip)
})

test_that("Validate Count Table returns a valid count table", {
  count_table <- read.delim(test_path("extdata", "amazon.count_table"))
  validated_count_table <- validate_count_table(count_table)

  expect_true(ncol(validated_count_table) == 3)
  expect_true("data.frame" %in% class(validated_count_table))
  expect_true("no_group" %in% colnames(validated_count_table))
  expect_false("no-group" %in% colnames(count_table))
})

test_that("Validate Count table does not change valid count tables", {
  count_table <- read.delim(test_path("extdata", "amazon.count_table"))
  validated_count_table <- validate_count_table(count_table)
  valid_count_table <- validate_count_table(validated_count_table)
  expect_true(ncol(valid_count_table) == ncol(valid_count_table))
})

test_that("Create sparse matrix will create valid sparse matrices", {
  i_values <- as.integer(1:100)
  j_values <- as.integer(sample(1:100, 100, TRUE))
  x_values <- as.numeric(runif(100, 0, 1))
  s_matrix <- create_sparse_matrix(i_values, j_values, x_values)
  expect_equal(length(s_matrix@i), length(i_values))
  expect_equal(length(s_matrix@j), length(j_values))
  expect_equal(length(s_matrix@x), length(x_values))
})

test_that("Create sparse matrix will fail when given invalid data", {
  i_values <- as.integer(1:101)
  j_values <- as.integer(sample(1:100, 100, TRUE))
  x_values <- as.numeric(runif(100, 0, 1))
  expect_error(create_sparse_matrix(i_values, j_values, x_values))
})

test_that("Split Clusters to list will generate valid list", {
  cutoff <- 0.2
  count_table <- read_count(test_path("extdata", "amazon.count_table"))
  distance_data <- read_dist(test_path("extdata", "amazon_column.dist"),
                             count_table, cutoff, FALSE)
  df <- cluster(distance_data)
  list <- split_clusters_to_list(df)
  expect_true(all(names(list) %in% df$cluster$otu))
  expect_true(length(list) == nrow(df$cluster))

  # Test that it works on classic cluster methods
  df <- cluster(distance_data, "furthest")
  list <- split_clusters_to_list(df)
  expect_true(all(names(list) %in% df$cluster$otu))
  expect_true(length(list) == nrow(df$cluster))
})
