# Pull the mass dataset and the distance
# Extract name file
# Extract count file
# muMS2 -> cluster file
# Samples contains names:
# Prototype is complete
Opticluster <- function(sparse_matrix, cutoff, iterations)
{
  indexOneList <- sparse_matrix@i
  indexTwoList <- sparse_matrix@j
  valueList <- sparse_matrix@x
  matrix_length <- length(indexOneList)
  clustering_output_string <- Opticluster::MatrixToOpiMatrixCluster(indexOneList, indexTwoList, valueList, cutoff,
                                   matrix_length, matrix_length, iterations)
  return(t(read.table(text = clustering_output_string, sep = "\t", header = TRUE)))
}