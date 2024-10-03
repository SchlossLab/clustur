#' clustur Description
#'
#' @export
#' @param cutoff A cutoff value
#' @param count_table A table of names and the given abundance per group.
#' @param iterations The number of iterations
#' @param shuffle a boolean to determine whether or
#'  not you want to shuffle the data before you cluster
#' @param simularity_matrix are you using a simularity matrix or distance matrix
#' @param random_seed you can set your own random
#' seed for consistent results, if not it will be set to 123
#' @param ... Either your phylip file or column file path,
#'  or a sparse distance matrix
#' @description
#' You must specfiy the type of matrix you are inputting
#'  to cluster your object and we support three types:
#' the path to your phylip and column distance file, or a sparse matrix.
#'
#' @examples
#'  # Using a sparse matrix
#'  i_values <- as.integer(1:100)
#'  j_values <- as.integer(sample(1:100, 100, TRUE))
#'  x_values <- as.numeric(runif(100, 0, 1))
#'  s_matrix <- Matrix::spMatrix(nrow=max(i_values),
#'                               ncol=max(i_values),
#'                               i=i_values,
#'                               j=j_values,
#'                               x=x_values)
#'
#'  # Creating a count table using the sparse matrix
#'  count_table_sparse <- data.frame(sequence=as.character(i_values),
#'                                  total=rep(1,times=100))
#'
#'  cluster_results <- opti_cluster(cutoff=0.2,
#'                                  count_table = count_table_sparse,
#'                                  sparse_matrix=s_matrix)
#'
#'  # With a column file
#'  count_table <- read.delim(example_path("amazon1.count_table"))
#'  amazon_data_column <- opti_cluster(column_path=
#'                                     example_path("96_sq_column_amazon.dist"),
#'                                     count_table = count_table, cutoff = 0.2)
#'  # With a phylip file
#'  count_table <- read.delim(example_path("amazon1.count_table"))
#'  amazon_data_phylip <- opti_cluster(phylip_path=
#'                                     example_path("98_sq_phylip_amazon.dist"),
#'                                     count_table = count_table, cutoff = 0.2)
#'
#'
#'
#' @return A data.frame of the cluster and cluster metrics.
opti_cluster <- function(cutoff, count_table,
                         iterations = 100, shuffle = TRUE,
                         simularity_matrix = FALSE, random_seed = 123, ...) {
  count_table <- validate_count_table(count_table)
  list_params <- list(...)
  params <- names(list_params)
  cluster_dfs <- list()

  if ("phylip_path" %in% params &&
        "column_path" %in% params &&
        "sparse_matrix" %in% params) {
    stop("You cannot use all three input paramters at once.
    Use either phylip_path, column_path, or sparse_matrix.")
  }
  set.seed(random_seed)
  if ("sparse_matrix" %in% params) {
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
  } else if ("phylip_path" %in% params) {
    phylip_path <- list_params$phylip_path
    cluster_dfs <- OptiClusterPhylip(
      phylip_path,
      cutoff,
      count_table,
      iterations,
      shuffle,
      simularity_matrix
    )
  } else if ("column_path" %in% params) {
    column_path <- list_params$column_path
    cluster_dfs <- OptiClusterColumnDist(
      column_path,
      cutoff,
      count_table,
      iterations,
      shuffle,
      simularity_matrix
    )
  } else {
    stop("The parameters should include either a sparse_matrix,
    phylip_path, column_path")
  }
  cluster_dfs[[4]]$comma_count <- sapply(cluster_dfs[[4]]$bins, function(x) {
    ls <- gregexpr(",", x, fixed = TRUE)[[1]]
    if (ls[[1]] == -1) {
      return(0)
    } else {
      return(length(ls))
    }
  })
  cluster_dfs[[4]] <- cluster_dfs[[4]][order(cluster_dfs[[4]]$comma_count,
                                             decreasing = TRUE), ]
  cluster_dfs[[4]] <- cluster_dfs[[4]][, 1:3]
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
#' @param random_seed you can set your own random seed
#'  for consistent results, if not it will be set to 123
#' @param ... Either your phylip file or column file path,
#'  or a sparse distance matrix
#' @description
#' You must specfiy the type of matrix you are inputting
#'  to cluster your object and we support three types:
#' the path to your phylip and column distance file, or a sparse matrix.
#' @return A string of the given cluster.
#'
#' @examples
#'  # Using a sparse matrix
#'  i_values <- as.integer(1:100)
#'  j_values <- as.integer(sample(1:100, 100, TRUE))
#'  x_values <- as.numeric(runif(100, 0, 1))
#'  s_matrix <- Matrix::spMatrix(nrow=max(i_values),
#'                               ncol=max(i_values),
#'                               i=i_values,
#'                               j=j_values,
#'                               x=x_values)
#'
#'  # Creating a count table using the sparse matrix
#'  count_table_sparse <- data.frame(sequence=as.character(i_values),
#'                                  total=rep(1,times=100))
#'  # furthest method
#'  cluster_results <- cluster(cutoff=0.2, count_table = count_table_sparse,
#'                             sparse_matrix=s_matrix, method="furthest")
#'
#'  # With a phylip file and nearest methods
#'  count_table <- read.delim(example_path("amazon1.count_table"))
#'  amazon_data_phylip <- cluster(phylip_path=
#'                                example_path("98_sq_phylip_amazon.dist"),
#'               count_table = count_table, method="nearest", cutoff = 0.2)
#'
#'  # With a column file and average methods
#'  amazon_data_column <- cluster(column_path=
#'                                example_path("96_sq_column_amazon.dist"),
#'               count_table = count_table, method="average", cutoff = 0.2)
#'
#'  # Weighted method
#'  amazon_data_column <- cluster(column_path=
#'                                example_path("96_sq_column_amazon.dist"),
#'              count_table = count_table, method="weighted", cutoff = 0.2)
#'
#'
cluster <- function(cutoff, method, count_table,
                    simularity_matrix = FALSE, random_seed = 123, ...) {
  list_params <- list(...)
  params <- names(list_params)
  cluster_dfs <- list()
  if ("phylip_path" %in% params &&
        "column_path" %in% params &&
        "sparse_matrix" %in% params) {
    stop("You cannot use all three input paramters at once.
    Use either phylip_path, column_path, or sparse_matrix.")
  }
  set.seed(random_seed)
  if ("sparse_matrix" %in% params) {
    sparse_matrix <- list_params$sparse_matrix
    cluster_dfs <-  ClassicCluster(
      sparse_matrix@i, sparse_matrix@j,
      sparse_matrix@x, cutoff, method,
      validate_count_table(count_table),
      simularity_matrix
    )
  } else if ("phylip_path" %in% params) {
    phylip_path <- list_params$phylip_path
    cluster_dfs <- ClusterWithPhylip(
      phylip_path, cutoff, method,
      validate_count_table(count_table),
      simularity_matrix
    )
  } else if ("column_path" %in% params) {
    column_path <- list_params$column_path
    cluster_dfs <-  ClusterWithColumn(
      column_path, cutoff, method,
      validate_count_table(count_table),
      simularity_matrix
    )
  } else {
    stop("The parameters should include either a sparse_matrix,
    phylip_path, column_path")
  }
  cluster_dfs[[2]]$comma_count <- sapply(cluster_dfs[[2]]$bins, function(x) {
    ls <- gregexpr(",", x, fixed = TRUE)[[1]]
    if (ls[[1]] == -1) {
      return(0)
    } else {
      return(length(ls))
    }
  })
  cluster_dfs[[2]] <- cluster_dfs[[2]][order(cluster_dfs[[2]]$comma_count,
                                             decreasing = TRUE), ]
  cluster_dfs[[2]] <- cluster_dfs[[2]][, 1:3]

  return(list(
    abundance = cluster_dfs[[1]],
    cluster = cluster_dfs[[2]]
  ))
}


validate_count_table <- function(count_table_df) {
  if (ncol(count_table_df) > 2) {
    return(count_table_df)
  }
  totals <- count_table_df$total
  count_table_df <- cbind(count_table_df, totals)
  names(count_table_df)[3] <- "no_group"
  count_table_df[[1]] <- as.character(count_table_df[[1]])
  return(count_table_df)
}


#' Example Path
#'
#' @export
#' This function was created as a helper function to generate file paths to our
#'  internal data. You are able to access this function if you
#'  want to follow along with the example.
#' @param file The data of the path you are looking to find.
#' @examples
#' # This will return the path to our example file
#' example_path("98_sq_phylip_amazon.dist")
#'
#' @return the path inside of the package of the file.
example_path <- function(file = NULL) {
  path <- ""
  if (is.null(file)) {
    path <- dir(system.file("extdata", package = "clustur"))
  } else {
    path <- system.file("extdata", file, package = "clustur", mustWork = TRUE)
  }
  return(path)
}


#' Convert Distance File To Sparse
#'
#' @export
#' This function will take your phylip or column distance file and convert it to a sparse matrix for easier use
#' in clustur
#' @param count_table Your count table based on the data.
#' @param distance_file_path The path to the data file.
#' @param type_of_file What type of file you are passing to the function, has to be "phylip", or "column".
#' @examples
#' # This will return the path to our example file
#' column_path <- example_path("updated_column.dist")
#' phylip_path <- example_path("updated_phylip_1.txt")
#' count_table <- readRDS(example_path("count_table.RDS"))
#' 
#' phylip_sparse <- convert_distance_file_to_sparse(count_table, phylip_path, "phylip")
#' column_sparse <- convert_distance_file_to_sparse(count_table, column_path, "column")
#' @return a sparse matrix based off the data provided
convert_distance_file_to_sparse <- function(count_table, distance_file_path,  type_of_file) {
  if(type_of_file != "phylip" && type_of_file != "column")
  {
    stop("The type_of_file parameter has to be 'phylip' or 'column.'")
  }
  # Have to convert c++ indexes to r indexes, so we add one to the i and j list
  triplicate_list <- DistanceFileToSparseMatrix(count_table, distance_file_path, type_of_file)
  named_triplicates <- list(i_values = triplicate_list[[1]] + 1, j_values = triplicate_list [[2]] + 1, data = triplicate_list[[3]])
  max_size <- max(named_triplicates$i_values, named_triplicates$j_values)
  return(spMatrix(nrow = max_size, ncol = max_size,
   i = named_triplicates$i_values, j = named_triplicates$j_values, x = named_triplicates$data))
}
