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
#' @param cutoff A cutoff value
#' @param count_table A table of names and the given abundance per group.
#' @param iterations The number of iterations
#' @param shuffle a boolean to determine whether or
#'  not you want to shuffle the data before you cluster
#' @param simularity_matrix are you using a simularity matrix or distance matrix
#' @param ... Either your phylip file or column file path, or a sparse distance matrix
#' @description
#' You must specfiy the type of matrix you are inputting to cluster your object and we support three types:
#' the path to your phylip and column distance file, or a sparse matrix.

#' @return A data.frame of the cluster and cluster metrics.
opti_cluster <- function(cutoff, count_table,
                         iterations = 100, shuffle = TRUE,
                         simularity_matrix = FALSE, ...) {
  count_table <- validate_count_table(count_table)
  list_params <- list(...)
  params <- names(list_params)
  cluster_dfs <- list()
  
  if("phylip_path" %in% params && 
    "column_path" %in% params &&
    "sparse_matrix" %in% params){
    stop("You cannot use all three input paramters at once.
    Use either phylip_path, column_path, or sparse_matrix.")
  }
  
  if("sparse_matrix" %in% params)
  {
    sparse_matrix <- list_params$sparse_matrix
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
  }
  else if("phylip_path" %in% params) {
    phylip_path <- list_params$phylip_path
    cluster_dfs <- OptiClusterPhylip(
      phylip_path,
      cutoff,
      count_table,
      iterations,
      shuffle,
      simularity_matrix
    )
  }
  else if("column_path" %in% params) {
    column_path <- list_params$column_path
    cluster_dfs <- OptiClusterColumnDist(
      column_path,
      cutoff,
      count_table,
      iterations,
      shuffle,
      simularity_matrix
    )
  }
  else {
    stop("The parameters should include either a sparse_matrix,
    phylip_path, column_path")
  }
  cluster_dfs[[4]]$comma_count <- sapply(cluster_dfs[[4]]$bins, function(x){
    ls <- gregexpr(",", x, fixed=TRUE)[[1]]
    if(ls[[1]] == -1){
      return(0)
    }
    else{
      return(length(ls))
    }
  })
  cluster_dfs[[4]] <- cluster_dfs[[4]][order(cluster_dfs[[4]]$comma_count, decreasing = T), ]
  cluster_dfs[[4]] <- cluster_dfs[[4]][,1:3]
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
#' @param cutoff A cutoff value.
#' @param method The type of cluster you wish to conduct;
#'  furthest, nearest, average, weighted.
#' @param count_table A table of names and the given abundance per group.
#' @param simularity_matrix are you using a simularity matrix or distance matrix
#' @param ... Either your phylip file or column file path, or a sparse distance matrix
#' @description
#' You must specfiy the type of matrix you are inputting to cluster your object and we support three types:
#' the path to your phylip and column distance file, or a sparse matrix.
#' @return A string of the given cluster.
cluster <- function(cutoff, method,
                    count_table, simularity_matrix = FALSE, ...) {
  list_params <- list(...)
  params <- names(list_params)
  cluster_dfs <- list()
  if("phylip_path" %in% params && 
    "column_path" %in% params &&
    "sparse_matrix" %in% params){
    stop("You cannot use all three input paramters at once.
    Use either phylip_path, column_path, or sparse_matrix.")
  }
  
  if("sparse_matrix" %in% params)
  {
    sparse_matrix <- list_params$sparse_matrix
    cluster_dfs <-  ClassicCluster(
      sparse_matrix@i, sparse_matrix@j,
      sparse_matrix@x, cutoff, method,
      validate_count_table(count_table),
      simularity_matrix
    )
  }
  else if("phylip_path" %in% params) {
    phylip_path <- list_params$phylip_path
    cluster_dfs <- ClusterWithPhylip(
      phylip_path, cutoff, method,
      validate_count_table(count_table),
      simularity_matrix
    )
  }
  else if("column_path" %in% params) {
    column_path <- list_params$column_path
    cluster_dfs <-  ClusterWithColumn(
      column_path, cutoff, method,
      validate_count_table(count_table),
      simularity_matrix
    )
  }
  else {
    stop("The parameters should include either a sparse_matrix,
    phylip_path, column_path")
  }
  
  cluster_dfs[[2]]$comma_count <- sapply(cluster_dfs[[2]]$bins, function(x){
    ls <- gregexpr(",", x, fixed=TRUE)[[1]]
    if(ls[[1]] == -1){
      return(0)
    }
    else{
      return(length(ls))
    }
  })
  cluster_dfs[[2]] <- cluster_dfs[[2]][order(cluster_dfs[[2]]$comma_count, decreasing = T), ]
  cluster_dfs[[2]] <- cluster_dfs[[2]][,1:3]

  return(list(
    abundance = cluster_dfs[[1]],
    cluster = cluster_dfs[[2]]
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
  order <- 1:nrow(count_table_df)
  count_table_df <- cbind(count_table_df, totals)
  count_table_df <- cbind(count_table_df, order)
  names(count_table_df)[3] <- "no_group"
  count_table_df[[1]] <- as.character(count_table_df[[1]])
  return(count_table_df)
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


example <- function(file = NULL) {
  path <- ""
  if (is.null(file)) {
    path <- dir(system.file("extdata", package = "Opticluster"))
  } else {
    path <- system.file("extdata", file, package = "Opticluster", mustWork = TRUE)
  }
  return(path)
}