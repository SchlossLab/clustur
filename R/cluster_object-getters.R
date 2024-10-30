#' Get Cutoff
#'
#' Returns the distance cutoff of the cluster object from the `cluster()`
#' function

#' @param cluster_data The output from the `cluster()` function.
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.full.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"),
#'                            count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opticlust")
#' cutoff <- get_cutoff(df_clusters)
#' @return the cutoff value as a `dbl`
#' @export

get_cutoff <- function(cluster_data) {
  if (!("list" %in% class(cluster_data)) ||
        is.null(cluster_data$label)) {
    stop("Ensure `cluster_data` is the object generated from the
    `cluster` function.")
  }
  return(cluster_data$label)
}


#' Get Clusters
#'
#' @export
#' @description
#' GetClusters returns a `data.frame` of the generated
#' clusters from the `cluster()` function.
#' @param cluster_data The output from the `cluster()` function.
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.full.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"),
#'                            count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opticlust")
#' clusters <- get_bins(df_clusters)
#'
#' @return the created cluster `data.frame`.
get_bins <- function(cluster_data) {
  if (!("list" %in% class(cluster_data)) ||
        is.null(cluster_data$cluster)) {
    stop("Ensure `cluster_data` is the object generated from the
    `cluster` function.")
  }
  return(cluster_data$cluster)
}


#' Get Shared
#'
#' @export
#' @description
#' GetShared retuns the generated abundance `data.frame`
#' from the `cluster()` function
#' @param cluster_data The output from the `cluster()` function.
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.full.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"),
#'                            count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opticlust")
#' shared <- get_abundance(df_clusters)
#'
#' @return a shared data.frame
get_abundance <- function(cluster_data) {
  if (!("list" %in% class(cluster_data)) ||
        is.null(cluster_data$abundance)) {
    stop("Ensure `cluster_data` is the object generated from the
    `cluster` function.")
  }
  return(cluster_data$abundance)
}


#' Get Metrics
#'
#' @export
#' @description
#' GetMetrics returns the generated metrics `data.frame` from the
#' `cluster()` function.
#' @param cluster_data The output from the `cluster()` function.
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.full.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"),
#'                            count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opticlust")
#' list_of_metrics <- get_metrics(df_clusters)
#'
#' @return a list of metric data.frames
get_metrics <- function(cluster_data) {
  if (!("list" %in% class(cluster_data)) ||
        is.null(cluster_data$cluster_metrics)) {
    stop("Can only use the get_metrics function when an 
    object is clustered using the opticluster method")
  }
  return(list(metrics = cluster_data$cluster_metrics,
              iteration_metrics = cluster_data$iteration_metrics))
}
