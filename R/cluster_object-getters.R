
#' Get Cutoff
#'
#' @export
#' @description
#' This function returns the cutoff of the cluster object
#' @param cluster_data The cluster method output
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"), count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opti")
#' cutoff <- get_cutoff(df_clusters)
#'
#' @return the cutoff value
get_cutoff <- function(cluster_data) {
  return(as.numeric(cluster_data$cluster$label[[1]]))
}


#' Get Clusters
#'
#' @export
#' @description
#' This function returns the cluster data.frame of the cluster object
#' @param cluster_data The cluster method output
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"), count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opti")
#' clusters <- get_clusters(df_clusters)
#'
#' @return the cluster data.frame
get_clusters <- function(cluster_data) {
  return(cluster_data$cluster)
}


#' Get Shared
#'
#' @export
#' @description
#' This function returns the shared data.frame of the cluster object
#' @param cluster_data The cluster method output
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"), count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opti")
#' shared <- get_shared(df_clusters)
#'
#' @return a shared data.frame
get_shared <- function(cluster_data) {
  return(cluster_data$abundance)
}


#' Get Metrics
#'
#' @export
#' @description
#' This function returns the all of the metrics from the cluster object
#' @param cluster_data The cluster method output
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"), count_table, cutoff, FALSE)
#' df_clusters <- cluster(distance_data, method = "opti")
#' list_of_metrics <- get_metrics(df_clusters)
#'
#' @return a list of metric data.frames
get_metrics <- function(cluster_data) {
  if(is.null(cluster_data$cluster_metrics)) {
    stop("Can only use the get_metrics function when an 
    object is clustered using the opticluster method")
  }
  return(list(metrics = cluster_data$cluster_metrics, 
              other_metrics = cluster_data$other_cluster_metrics))
}

