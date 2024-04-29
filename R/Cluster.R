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
#' @param param1 A Sparse Matrix.
#' @param param2 A cutoff value
#' @param param3 The number of iterations
#' @return A data.frame of the clusters.
OptiCluster <- function(sparse_matrix, cutoff, iterations)
{
  indexOneList <- sparse_matrix@i
  indexTwoList <- sparse_matrix@j
  valueList <- sparse_matrix@x
  matrix_length <- length(indexOneList)
  clustering_output_string <- MatrixToOpiMatrixCluster(indexOneList, indexTwoList, valueList, cutoff,
                                   matrix_length, matrix_length, iterations)
  df <- t(read.table(text = clustering_output_string, sep = "\t", header = TRUE))
  colnames(df)[1] <- "cluster"
  return(data.frame(df))
}