#' Split Clusters to List
#'
#' @export
#' @description
#' `split_clusters_to_list()` will extract clusters from the
#' cluster generated `data.frame`. It will then turn those clusters
#' into a list. This allows users to more easily visualize their data.
#' @param cluster The output generated from the `cluster()` function.
#' @examples
#'
#'  cutoff <- 0.2
#'  count_table <- read_count(example_path("amazon.full.count_table"))
#'  distance_data <- read_dist(example_path("amazon_column.dist"),
#'                             count_table, cutoff, FALSE)
#'  cluster_results <- cluster(distance_data, method = "opticlust")
#'
#'  cluster_list <- split_clusters_to_list(cluster_results)
#' @return a named `list` of clusters.
split_clusters_to_list <- function(cluster) {
  ls <- c()
  cluster_df <- get_clusters(cluster)
  size <- nrow(cluster_df)
  for (i in 1:size) {
    ls[[cluster_df$otu[[i]]]] <-
      as.list(el(strsplit(cluster_df$sequences[[i]], ",")))
  }
  return(ls)
}
