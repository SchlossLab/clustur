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
#' @param count_table A table of names and the given abundance per group.
#' @param iterations The number of iterations
#' @param shuffle a boolean to determine whether or
#'  not you want to shuffle the data before you cluster
#' @param simularity_matrix are you using a simularity matrix or distance matrix
#' @return A data.frame of the cluster and cluster metrics.
opti_cluster <- function(sparse_matrix, cutoff, count_table,
                         iterations = 100, shuffle = TRUE,
                         simularity_matrix = FALSE) {
  count_table <- validate_count_table(count_table)
  cluster_dfs <- MatrixToOpiMatrixCluster(
    sparse_matrix@i,
    sparse_matrix@j,
    sparse_matrix@x,
    cutoff,
    count_table,
    iterations,
    shuffle,
    simularity_matrix
  )
  opticluster_data <- list(
    abundance = cluster_dfs[[1]],
    cluster = cluster_dfs[[4]],
    cluster_metrics = cluster_dfs[[3]],
    other_cluster_metrics = cluster_dfs[[2]]
  )

  return(opticluster_data)
}

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
#' @param phylip_path The path of your phylip file
#' @param cutoff A cutoff value
#' @param count_table A table of names and the given abundance per group.
#' @param iterations The number of iterations
#' @param shuffle a boolean to determine whether or
#'  not you want to shuffle the data before you cluster
#' @param simularity_matrix are you using a simularity matrix or distance matrix
#' @return A data.frame of the cluster and cluster metrics.
opti_cluster_phylip <- function(phylip_path, cutoff, count_table,
  iterations = 100, shuffle = TRUE,
  simularity_matrix = FALSE) {
count_table <- validate_count_table(count_table)
cluster_dfs <- OptiClusterPhylip(
phylip_path,
cutoff,
count_table,
iterations,
shuffle,
simularity_matrix
)
opticluster_data <- list(
abundance = cluster_dfs[[1]],
cluster = cluster_dfs[[4]],
cluster_metrics = cluster_dfs[[3]],
other_cluster_metrics = cluster_dfs[[2]]
)

return(opticluster_data)
}

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
#' @param column_path The path of your phylip file
#' @param cutoff A cutoff value
#' @param count_table A table of names and the given abundance per group.
#' @param iterations The number of iterations
#' @param shuffle a boolean to determine whether or
#'  not you want to shuffle the data before you cluster
#' @param simularity_matrix are you using a simularity matrix or distance matrix
#' @return A data.frame of the cluster and cluster metrics.
opti_cluster_column <- function(column_path, cutoff, count_table,
  iterations = 100, shuffle = TRUE,
  simularity_matrix = FALSE) {
count_table <- validate_count_table(count_table)
cluster_dfs <- OptiClusterColumnDist(
column_path,
cutoff,
count_table,
iterations,
shuffle,
simularity_matrix
)
opticluster_data <- list(
abundance = cluster_dfs[[1]],
cluster = cluster_dfs[[4]],
cluster_metrics = cluster_dfs[[3]],
other_cluster_metrics = cluster_dfs[[2]]
)

return(opticluster_data)
}




#' Cluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param sparse_matrix A Sparse Matrix.
#' @param cutoff A cutoff value.
#' @param method The type of cluster you wish to conduct;
#'  furthest, nearest, average, weighted.
#' @param count_table A table of names and the given abundance per group.
#' @param simularity_matrix are you using a simularity matrix or
#'  distance matrix.
#' @return A string of the given cluster.
cluster <- function(sparse_matrix, cutoff, method,
                    count_table, simularity_matrix = FALSE) {
  df <- ClassicCluster(
    sparse_matrix@i, sparse_matrix@j,
    sparse_matrix@x, cutoff, method,
    validate_count_table(count_table),
    simularity_matrix
  )
  return(list(
    abundance = df[[1]],
    cluster = df[[2]]
  ))
}


#' Cluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param phylip_path A phylip file path.
#' @param cutoff A cutoff value.
#' @param method The type of cluster you wish to conduct;
#'  furthest, nearest, average, weighted.
#' @param count_table A table of names and the given abundance per group.
#' @param simularity_matrix are you using a simularity matrix or
#'  distance matrix.
#' @return A string of the given cluster.
clusterPhylip <- function(phylip_path, cutoff, method,
  count_table, simularity_matrix = FALSE) {
df <- ClusterWithPhylip(
  phylip_path, 
  cutoff, method,
  validate_count_table(count_table),
  simularity_matrix
)
return(list(
abundance = df[[1]],
cluster = df[[2]]
))
}


#' Cluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param column_matrix_path A phylip file path.
#' @param cutoff A cutoff value.
#' @param method The type of cluster you wish to conduct;
#'  furthest, nearest, average, weighted.
#' @param count_table A table of names and the given abundance per group.
#' @param simularity_matrix are you using a simularity matrix or
#'  distance matrix.
#' @return A string of the given cluster.
clusterColumn <- function(column_matrix_path, cutoff, method,
  count_table, simularity_matrix = FALSE) {
df <- ClusterWithColumn(
  column_matrix_path, 
  cutoff, method,
  validate_count_table(count_table),
  simularity_matrix
)
return(list(
abundance = df[[1]],
cluster = df[[2]]
)) 
}


#' Opticluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param count_table_df The count table,
#' which contains all your abundance information for each sequence.
#' @return Validated count table.
validate_count_table <- function(count_table_df) {
  if (ncol(count_table_df) > 2) {
    return(count_table_df)
  }
  totals <- count_table_df$total
  count_table_df <- cbind(count_table_df, totals)
  names(count_table_df)[3] <- "no_group"
  count_table_df[[1]] <- as.character(count_table_df[[1]])
}

clust.abc <- function(x, ...)
{
    print("abc")
}

clust.bc <- function(x, ...)
{
  print("bcd")
}
clust <- function(x, ...)
{
  print(x)
}

