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
  clustering_output_string <- MatrixToOpiMatrixCluster(index_one_list, index_two_list, value_list, cutoff,
                                                       iterations, shuffle)
  df <- t(read.table(text = clustering_output_string,
                     sep = "\t", header = TRUE))
  df <- data.frame(df[-1, ])
  colnames(df)[1] <- "cluster"
  return(df)
}
