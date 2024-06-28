# Pull the mass dataset and the distance
# Extract name file
# Extract count file
# muMS2 -> cluster file
# Samples contains names:
# Prototype is complete
#' Opticluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param sparse_matrix A Sparse Matrix.
#' @param cutoff A cutoff value
#' @param iterations The number of iterations
#' @param shuffle a boolean to determine whether or not you want to shuffle the data before you cluster
#' @return A data.frame of the clusters.
opti_cluster <- function(sparse_matrix, cutoff, iterations, shuffle = TRUE) {
  index_one_list <- sparse_matrix@i
  index_two_list <- sparse_matrix@j
  value_list <- sparse_matrix@x
  clustering_output_string_list <- MatrixToOpiMatrixCluster(index_one_list, index_two_list, value_list, cutoff,
                                                       iterations, shuffle)
  clustering_output_string <- clustering_output_string_list[1]
  clustering_metric <- clustering_output_string_list[2]
  clustering_metric_2 <- clustering_output_string_list[3]
  df_cluster_metrics <- (read.table(text = clustering_metric,
                     sep = "\t", header = TRUE))
  df_other_cluster_metrics <-  (read.table(text = clustering_metric_2,
                     sep = "\t", header = TRUE))

  df_cluster <- t(read.table(text = clustering_output_string,
                     sep = "\t", header = TRUE))
  df_cluster <- data.frame(df_cluster[-1, ])

  colnames(df_cluster)[1] <- "cluster"

  opticluster_data <- list(cluster = df_cluster,
                           cluster_metrics = df_cluster_metrics,
                           other_cluster_metrics = df_other_cluster_metrics)

  return(opticluster_data)
}

create_phylip_file <- function(sparse_matrix, cutoff, file_path)
{
  return (MatrixDataToPhylipFormat(sparse_matrix@i, sparse_matrix@j,
                           sparse_matrix@x, cutoff, file_path))
}
cluster <- function(sparse_matrix, cutoff, method)
{
  return (ClassicCluster(sparse_matrix@i, sparse_matrix@j,
                           sparse_matrix@x, cutoff, method))
}
function_1 <- function(){
  devtools::load_all()
    matrix <- readRDS(test_path("extdata","matrix_data.RDS"))
    # dat_string <- create_phylip_file(matrix, 0.2,
    #                "/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/matrix_output1000.txt")
  clust <- cluster(matrix, 0.2, "furthest")
  write(clust, "tests/testthat/extdata/Furthest_full.txt")

}