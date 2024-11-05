#' Get Count Table
#'
#' @export
#' @description
#' This function returns the count table that was used
#' to generate the distance object.
#' @param distance_object The output from the `read.dist()` function.
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.full.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"),
#'                            count_table, cutoff, FALSE)
#' count_table <- get_count_table(distance_data)
#'
#' @return a count_table `data.frame`.
get_count_table <- function(distance_object) {
  if (!inherits(distance_object, "externalptr"))
    stop("The distance object must be the object generated when
         calling the `read_dist` function.")
  return(GetCountTable(distance_object))
}

#' Get Distance Data Frame
#'
#' @export
#' @description
#' This function will generate a `data.frame` that contains the distances
#' of all the indexes.
#' @param distance_object The output from the `read.dist()` function.
#' @examples
#' cutoff <- 0.2
#' count_table <- read_count(example_path("amazon.full.count_table"))
#' distance_data <- read_dist(example_path("amazon_column.dist"),
#'                                         count_table, cutoff, FALSE)
#' count_table <- get_count_table(distance_data)
#'
#' @return a distance `data.frame`.
get_distance_df <- function(distance_object) {
  if (!inherits(distance_object, "externalptr"))
    stop("The distance object must be the object generated when
         calling the `read_dist` function.")
  df <- GetDistanceDataFrame(distance_object)
  return(df)
}


# list_to_remove <- c()
# for(i in 1:nrow(dist_df_sparse)) {
#   first_name <- dist_df_sparse[i,]$FirstName
#   second_name <- dist_df_sparse[i,]$SecondName
#   my_index <- which(dist_df_sparse$FirstName == first_name & 
#     dist_df_sparse$SecondName == second_name, arr.ind = TRUE)

#   index <- which(dist_df_sparse$FirstName == second_name & 
#     dist_df_sparse$SecondName == first_name, arr.ind = TRUE)
#   if(length(index) > 0 && !(my_index %in% list_to_remove)) {
#     list_to_remove <- c(list_to_remove, index)
#   }
# }
# dist_df_sparse[-c(list_to_remove),]
