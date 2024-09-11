test_that("Clustering returns proper results", {
  expected_df <- readRDS(test_path("extdata","abundanceResult.RDS"))
  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  df <- Opticluster::opti_cluster(sparse_matrix, 0.2, count_table)
  # trueList <- do.call(paste0, df$abundance) %in% do.call(paste0, expected_df)
  expect_equal(class(df$cluster), "data.frame")
  expect_equal(class(df$cluster_metrics), "data.frame")
  expect_equal(class(df$other_cluster_metrics), "data.frame")
  expect_equal(class(df$abundance), "data.frame")
  # expect_true(all(trueList == TRUE))
})

#Reader function
# label, ids, bins
# opt to remove label

# For cluster
# Sequence, Abundance, OTU
test_that("Normal Cluster is able to properly cluster data",
{

  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))

  cluster_furthest <- cluster(sparse_matrix, 0.2, "furthest", count_table, FALSE)
  cluster_average <- cluster(sparse_matrix, 0.2, "average", count_table, FALSE) ## Might actually be a rounding issue...????
  cluster_weighted <- cluster(sparse_matrix, 0.2, "weighted", count_table, FALSE) # TODO Fix Label
  cluster_nearest <- cluster(sparse_matrix, 0.2, "nearest", count_table, FALSE)

  expect_length(cluster_furthest, 2)
  expect_length(cluster_average, 2)
  expect_length(cluster_weighted, 2)
  expect_length(cluster_nearest, 2)


  expect_true(any(class(cluster_furthest$abundance) == "data.frame"))
  expect_true(any(class(cluster_average$abundance) == "data.frame"))
  expect_true(any(class(cluster_weighted$abundance) == "data.frame"))
  expect_true(any(class(cluster_nearest$abundance) == "data.frame"))
})

# tidy_data <- reshape2::melt(shared_file, id.vars=c("label" ), variable.name = "OTU" ,value.name = "abundance")
# cluster(phylip=/Users/grejoh/Documents/mothur/mothur/Clustur/updated_phylip.txt, count=/Users/grejoh/Documents/mothur/mothur/Clustur/count_table.count, cutoff=0.2)
