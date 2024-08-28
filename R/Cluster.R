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
  clustering_output_string_list <- MatrixToOpiMatrixCluster(index_one_list, index_two_list, value_list, cutoff,
                                                       iterations, shuffle)
  clustering_output_string <- clustering_output_string_list[1]
  clustering_metric <- clustering_output_string_list[2]
  clustering_metric_2 <- clustering_output_string_list[3]
  df_cluster_metrics <- (read.table(text = clustering_metric,
                     sep = "\t", header = TRUE))
  df_other_cluster_metrics <-  (read.table(text = clustering_metric_2,
                     sep = "\t", header = TRUE))

  df_cluster <- t(read.table(text = clustering_output_string,
                     sep = "\t", header = TRUE))
  df_cluster <- data.frame(df_cluster[-1, ])

  colnames(df_cluster)[1] <- "cluster"

  opticluster_data <- list(cluster = df_cluster,
                           cluster_metrics = df_cluster_metrics,
                           other_cluster_metrics = df_other_cluster_metrics)

  return(opticluster_data)
}



#' Opticluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param sparse_matrix A Sparse Matrix.
#' @param cutoff A cutoff value
#' @param method The type of cluster you wish to conduct. There are four different types:
#' furthest, nearest, average, weighted.
#' @return A string of the given cluster.
cluster <- function(sparse_matrix, cutoff, method)
{

  return (ClassicCluster(sparse_matrix@i, sparse_matrix@j,
                           sparse_matrix@x, cutoff, method, sampled_peak_table))
}

# df_read_table <- (read.table(text = cluster_furthest,
#   sep = "\t", header = TRUE, na.strings = "#NA", fill = TRUE))


# May have to create a formater package in the utils
progenesis_formatter <- function(peak_table) {
  peak_table <- readr::read_csv(peak_table,
    skip = 2,
    show_col_types = FALSE
  )
  names(peak_table)[2] <- "mz"
  names(peak_table)[3] <- "rt"
  return(peak_table)
}

# Simple way, but is hard coded. Does not account for other columns
create_count_table <- function(x)
{
  df <- progenesis_formatter("tests/testthat/extdata/102623_peaktable_coculture_simple.csv")
  df$total <- rowSums(df[which(!(names(df) == "Compound"))])
  df <- df[,which(!(names(df) %in% c("mz", "rt")))]
  df <- data.frame("Compound" = df$Compound, "total" = df$total, df[,which(!(names(df) %in% c("total", "Compound")))])
  meta_data <- read.csv("tests/testthat/extdata/102623_metadata_simple.csv")
  samples <- which((colnames(df) %in% grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)))
  other_columns <- names(df)[-samples]
  tidy_df <- melt(df, id = "Compound", variable.name = "sample", value.name = "intensity")
  groups <- names(df)[samples]
  # grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)
}

# In our case, we want to cluster only ms2 spectra, but we will force the user to give us a peak table that only has ms2 spectra.
# Other it clusters as needed
create_count_table_no_total <- function(x)
{
  df <- progenesis_formatter("tests\\testthat\\extdata\\102623_peaktable_coculture_simple.csv")
  df <- df[,which(!(names(df) %in% c("mz", "rt")))]
  meta_data <- read.csv("tests/testthat/extdata/102623_metadata_simple.csv")
  samples <- which((colnames(df) %in% grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)))
  other_columns <- names(df)[-samples]
  tidy_df <- melt(df, id = "Compound", variable.name = "sample", value.name = "intensity")
  groups <- names(df)[samples]
  # grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)
}

create_fake_sparse_matrix <- function(peak_table)
{
  peak_table <- df
  count <- nrow(peak_table)
  number_of_ions <- ceiling(runif(n=1, min=count/2, max=count))
  number_of_ions <- 1046
  sampled_peak_table <- df[sample(nrow(df), number_of_ions), ] # Ions with ms2 peaks or the count_table
  I_values <- sort(sampled_peak_table$Compound, decreasing = FALSE)
  Y_values <- sampled_peak_table$Compound
  values <- runif(number_of_ions, .01, .99)
  T2 <- new("dgTMatrix",
          i = as.integer(I_values),
          j = as.integer(Y_values), x=as.double(values), Dim=1325:1326)
  WritePhylipFile(I_values, Y_values, values, 0.25, "Phylip_file.txt")
}
# T2 <- new("dgTMatrix",
#           i = as.integer(c(1,1,0,3,3)),
#           j = as.integer(c(2,2,4,0,0)), x=10*1:5, Dim=1325:1326)
#   T2@x <- I_values

# sparse <- Matrix::sparseMatrix( i = I_values, j=Y_values, x = values)
# WritePhylipFile(I_values, Y_values, values, 0.25, "Phylip_file.txt")
# cluster(count=/Users/grejoh/Documents/OptiClusterPackage/Opticluster/count_table_sample.count_table, phylip=/Users/grejoh/Documents/OptiClusterPackage/Opticluster/Phylip_file.txt)
# # write
#  Need to figure out how to create a sparse distance matrix,
# Whether it be them creating a sparse distance matrix and inputting it as a phylip file
# Or using the scoring package to score things inside of this package.
# I think I will have them input a sparse distance file.
# Return shared file in tidy format and compare the result to mothur

# MOTHUR COMMAND
# cluster(count=/Users/grejoh/Documents/OptiClusterPackage/Opticluster/count_table_sample.count_table, phylip=/Users/grejoh/Documents/OptiClusterPackage/Opticluster/Phylip_file.txt)
# clust_fur <- ClassicCluster(I_values, Y_values,
#   values, 0.3, "furthest", sampled_peak_table)

# Create Count Table
# write.table(sampled_peak_table, "count_table_sample.count_table", quote = F, col.names = FALSE)
