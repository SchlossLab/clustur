test_that("Clustering returns proper results", {
  expected_df <- readRDS(test_path("extdata", "abundanceResult.RDS"))
  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  cutoff <- 0.2
  df <- Opticluster::opti_cluster(sparse_matrix, cutoff, count_table)
  expect_equal(class(df$cluster), "data.frame")
  expect_equal(class(df$cluster_metrics), "data.frame")
  expect_equal(class(df$other_cluster_metrics), "data.frame")
  expect_equal(class(df$abundance), "data.frame")
})

# Reader function
# label, ids, bins
# opt to remove label

# For cluster
# Sequence, Abundance, OTU
test_that("Normal Cluster is able to properly cluster data", {
  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  cluster_furthest <- cluster(
    sparse_matrix, 0.2,
    "furthest", count_table, FALSE
  )
  cluster_average <- cluster(
    sparse_matrix, 0.2,
    "average", count_table, FALSE
  )
  cluster_weighted <- cluster(
    sparse_matrix, 0.2,
    "weighted", count_table, FALSE
  )
  cluster_nearest <- cluster(
    sparse_matrix, 0.2,
    "nearest", count_table, FALSE
  )

  expect_length(cluster_furthest, 2)
  expect_length(cluster_average, 2)
  expect_length(cluster_weighted, 2)
  expect_length(cluster_nearest, 2)


  expect_true(any(class(cluster_furthest$abundance) == "data.frame"))
  expect_true(any(class(cluster_average$abundance) == "data.frame"))
  expect_true(any(class(cluster_weighted$abundance) == "data.frame"))
  expect_true(any(class(cluster_nearest$abundance) == "data.frame"))
})

test_that("Normal cluster works via phylip file", {

  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  phylip_path <- test_path("extdata", "updated_phylip_1.txt")
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  cluster_phylip_ff <- clusterPhylip(phylip_path, 
    0.2, "furthest",
    count_table, 
    FALSE)
  cluster_furthest <- cluster(
    sparse_matrix, 0.2,
    "furthest", count_table, FALSE
  )
  expect_true(all(cluster_furthest$cluster == cluster_phylip_ff$cluster))
  expect_true(all(cluster_furthest$abundance == cluster_phylip_ff$abundance))
})

test_that("Opticluster cluster works via phylip file", {
  df_column <- clusterColumn("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/column.txt", 0.2, "furthest", count_table, F)
  ClusterWithColumn("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/column.txt", 0.2, "furthest", count_table, F)
  expected_df <- readRDS(test_path("extdata", "abundanceResult.RDS"))
  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  phylip_path <- test_path("extdata", "updated_phylip_1.txt")
  cutoff <- 0.2
  df <- Opticluster::opti_cluster(sparse_matrix, cutoff, count_table)
  df2 <- opti_cluster_phylip(phylip_path, cutoff, count_table)
  expect_true(all(df$cluster == df2$cluster))
  expect_true(all(df$abundance == df2$abundance))
})
