test_that("Opticluster returns proper results", {
  expected_df <- readRDS(test_path("extdata", "abundanceResult.RDS"))
  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  cutoff <- 0.2
  df <- Opticluster::opti_cluster(sparse_matrix=sparse_matrix, cutoff, count_table)
  expect_equal(class(df$cluster), "data.frame")
  expect_equal(class(df$cluster_metrics), "data.frame")
  expect_equal(class(df$other_cluster_metrics), "data.frame")
  expect_equal(class(df$abundance), "data.frame")
})

test_that("Opticluster cluster works via phylip file", {

  expected_df <- readRDS(test_path("extdata", "abundanceResult.RDS"))
  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  phylip_path <- test_path("extdata", "updated_phylip_1.txt")
  cutoff <- 0.2
  df <- Opticluster::opti_cluster(sparse_matrix=sparse_matrix, cutoff, count_table)
  df2 <- Opticluster::opti_cluster(phylip_path=phylip_path, cutoff, count_table)
  expect_true(all(df$cluster == df2$cluster))
  expect_true(all(df$abundance == df2$abundance))
})

test_that("Normal cluster works via column file", {

  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  column_path <- test_path("extdata", "column.dist")
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))
  opticluster_column <- opti_cluster(
    column_path=column_path, 0.2,
    count_table, FALSE
  )
  opticluster_sparse <- opti_cluster(
    sparse_matrix=sparse_matrix, 0.2,
    count_table, FALSE
  )
  expect_true(all(opticluster_column$cluster$bins %in% opticluster_sparse$cluster$bins))
  expect_true(all(opticluster_column$abundance$abundance %in% opticluster_sparse$abundance$abundance))
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
    sparse_matrix=sparse_matrix, 0.2,
    "furthest", count_table, FALSE
  )
  cluster_average <- cluster(
    sparse_matrix=sparse_matrix, 0.2,
    "average", count_table, FALSE
  )
  cluster_weighted <- cluster(
    sparse_matrix=sparse_matrix, 0.2,
    "weighted", count_table, FALSE
  )
  cluster_nearest <- cluster(
    sparse_matrix=sparse_matrix, 0.2,
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
  cluster_phylip <- cluster(phylip_path=phylip_path, 
    0.2, "furthest",
    count_table, 
    FALSE)
  cluster_furthest <- cluster(
    sparse_matrix=sparse_matrix, 0.2,
    "furthest", count_table, FALSE
  )
  expect_true(all(cluster_furthest$cluster == cluster_phylip$cluster))
  expect_true(all(cluster_furthest$abundance == cluster_phylip$abundance))
})

test_that("Normal cluster works via column file", {

  sparse_matrix <- readRDS(test_path("extdata", "sparse_matrix.RDS"))
  column_path <- test_path("extdata", "column.dist")
  count_table <- readRDS(test_path("extdata", "count_table.RDS"))

  cluster_furthest_column<- cluster(
    column_path=column_path, 0.2, 
    "furthest", count_table, FALSE
  )
  expect_true(all(as.numeric(cluster_furthest_column$cluster$label) == 0))
})

test_that("Amazon Data from mothur clusters properly", {
  df_count <- readRDS(test_path("extdata", "amazon_count_table.RDS"))
  column_path <- test_path("extdata", "96_sq_column_amazon.dist")
  data <- opti_cluster(column_path=column_path, 
                        count_table = df_count, cutoff = 0.2)
  tidy_answer <- readRDS(test_path("extdata", "amazonResults.RDS"))
 
  # all(tidy_answer$clusters %in% 
  #   data$cluster$bins)
  ls <- as.list(el(strsplit(data$cluster$bins[[1]], ","))) 
  ls_answer <- as.list(el(strsplit(tidy_answer$cluster[[1]], ","))) 
  expect_true(all(ls %in% ls_answer))
  expect_true(length(tidy_answer$clusters) ==length(data$cluster$bins))
})
