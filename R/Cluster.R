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
#' @param count_table A table of abundances
#' @return A data.frame of the clusters.
opti_cluster <- function(sparse_matrix, cutoff, iterations = 100, shuffle = TRUE, count_table) {
  index_one_list <- sparse_matrix@i
  index_two_list <- sparse_matrix@j
  value_list <- sparse_matrix@x
  count_table <- validate_count_table(count_table)
  cluster_dfs <- MatrixToOpiMatrixCluster(index_one_list, index_two_list, value_list, cutoff,
                                                       count_table, iterations, shuffle)
  opticluster_data <- list(abundance = cluster_dfs[[1]],
                           cluster = cluster_dfs[[4]],
                           cluster_metrics = cluster_dfs[[3]],
                           other_cluster_metrics = cluster_dfs[[2]])

  return(opticluster_data)
}





#' Opticluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param sparse_matrix A Sparse Matrix.
#' @param cutoff A cutoff value
#' @param count_table A table of abundances
#' @param method The type of cluster you wish to conduct. There are four different types:
#' furthest, nearest, average, weighted.
#' @return A string of the given cluster.
cluster <- function(sparse_matrix, cutoff, method, simularity_matrix = FALSE, count_table)
{
  df <- ClassicCluster(sparse_matrix@i, sparse_matrix@j,
                           sparse_matrix@x, cutoff, method, simularity_matrix, validate_count_table(count_table))
  
  return (list(abundance = df[[1]],
                cluster = df[[2]]))
}


#' Opticluster Description
#'
#' Detailed description of the function.
#'
#' @export
#' @param count_table_df The count table, which contains all your abundance information for each sequence.
#' @return Validated count table.
validate_count_table <- function(count_table_df){
  if(ncol(count_table_df) > 2)
    return(count_table_df)
  totals <- count_table_df$total
  count_table_df <- cbind(count_table_df, totals)
  names(count_table_df)[3] <- "no_group"
  count_table_df[[1]] <- as.character(count_table_df[[1]])
  return(count_table_df)
}
# df_read_table <- (read.table(text = cluster_furthest,
#   sep = "\t", header = TRUE, na.strings = "#NA", fill = TRUE))


# # May have to create a formater package in the utils
# progenesis_formatter <- function(peak_table) {
#   peak_table <- readr::read_csv(peak_table,
#     skip = 2,
#     show_col_types = FALSE
#   )
#   names(peak_table)[2] <- "mz"
#   names(peak_table)[3] <- "rt"
#   return(peak_table)
# }

# # Simple way, but is hard coded. Does not account for other columns
# create_count_table <- function(x)
# {
#   df <- progenesis_formatter("tests/testthat/extdata/102623_peaktable_coculture_simple.csv")
#   df$total <- rowSums(df[which(!(names(df) == "Compound"))])
#   df <- df[,which(!(names(df) %in% c("mz", "rt")))]
#   df <- data.frame("Compound" = df$Compound, "total" = df$total, df[,which(!(names(df) %in% c("total", "Compound")))])
#   meta_data <- read.csv("tests/testthat/extdata/102623_metadata_simple.csv")
#   samples <- which((colnames(df) %in% grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)))
#   other_columns <- names(df)[-samples]
#   tidy_df <- melt(df, id = "Compound", variable.name = "sample", value.name = "intensity")
#   groups <- names(df)[samples]
#   # grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)
# }

# # In our case, we want to cluster only ms2 spectra, but we will force the user to give us a peak table that only has ms2 spectra.
# # Other it clusters as needed
# create_count_table_no_total <- function(x)
# {
#   df <- progenesis_formatter("tests\\testthat\\extdata\\102623_peaktable_coculture_simple.csv")
#   df <- df[,which(!(names(df) %in% c("mz", "rt")))]
#   meta_data <- read.csv("tests/testthat/extdata/102623_metadata_simple.csv")
#   samples <- which((colnames(df) %in% grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)))
#   other_columns <- names(df)[-samples]
#   tidy_df <- reshape(df, idvar = "Compound", direction="long")
#   groups <- names(df)[samples]
#   # grep(paste(meta_data$Sample_Code, collapse ="|"), names(df), value = TRUE)
# }

# create_fake_sparse_matrix <- function(peak_table)
# {
#   peak_table <- df
#   count <- nrow(peak_table)
#   number_of_ions <- ceiling(runif(n=1, min=count/2, max=count))
#   number_of_ions <- 1046
#   sampled_peak_table <- df[sample(nrow(df), number_of_ions), ] # Ions with ms2 peaks or the count_table
#   I_values <- sort(sampled_peak_table$Compound, decreasing = FALSE)
#   Y_values <- sampled_peak_table$Compound
#   values <- runif(number_of_ions, .01, .99)
#   max(Y_values)
#   T2 <- new("dgTMatrix",
#           i = as.integer(I_values),
#           j = as.integer(Y_values), x=as.double(values), Dim=1334:1335)
#   WritePhylipFile(I_values, Y_values, values, 0.25, "Phylip_file.txt")
# }



# # count_table_no_groups <- sampled_peak_table[, 1:2]
# # # count_table_no_groups <- validate_count_table(count_table_no_groups)
# T2 <- new("dgTMatrix",
#           i = as.integer(c(0,1,3,4,5)),
#           j = as.integer(c(1,0,3,4,5)),x=10*1:5, Dim=6:7)

# T2 <- new("dgTMatrix",
#           i = as.integer(c(1,1,0,3,3)),
#           j = as.integer(c(2,2,4,0,0)), x=10*1:5, Dim=4:5)


# T2@i <- as.integer(I_values)
# T2@j <- as.integer(Y_values)
# T2@x <- values

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
# # cluster(count=/Users/grejoh/Documents/OptiClusterPackage/Opticluster/count_table_sample.count_table, phylip=/Users/grejoh/Documents/OptiClusterPackage/Opticluster/Phylip_file.txt)
# clust_fur <- ClassicCluster(I_values, Y_values,
#    values, 0.3, "furthest", sampled_peak_table)
# sampled_peak_table$Compound <- as.character(sampled_peak_table$Compound)
# count_table_ms <- sampled_peak_table
# clust_fur <- cluster(T2, 0.3, "furthest", count_table_ms)
# clust_fur_no_groups <- cluster(T2, 0.3, "furthest", count_table_no_groups)

# TODO fix the opticluster return statement!!
# index_one_list <- T2@i
# index_two_list <- T2@j
# value_list <- T2@x
# count_table <- validate_count_table(count_table_no_groups)
# df <- MatrixToOpiMatrixCluster(index_one_list, index_two_list, value_list, 0.3,
#                                                      count_table, 1, TRUE)
# # # Create Count Table
# write.table(sampled_peak_table, "count_table_sample.count_table", quote = F, col.names = FALSE)
# tidy_shared <- readr::read_tsv("SharedFile_tidy.txt")



# table <- readRDS("test_table.RDS")
# T2 <- new("dgTMatrix",
#           i = as.integer(c(0,1,3,4,5)),
#           j = as.integer(c(1,0,3,4,5)),x=10*1:5, Dim=6:7)
# table <- test
# i_values  <- c()
# y_values <- c()
# data <-c()
# for(i in 2:ncol(table))
# {
#   i_values <- c(i_values, rep(i - 1, 10))
#   y_values <- c(y_values, 1:10)
#   data <- c(data, table[[i]])
# }

# T2@i <- as.integer(i_values)
# T2@j <- y_values
# T2@x <- data
#saveRDS(T2, "sparse_amazon.RDS")

# table$Names
 #sparse_amazon <- readRDS("sparse_amazon.RDS")
# saveRDS(table, "test_table.RDS")
# saveRDS(count <- data.frame(Names = table$Names, total = rep(1, 10)), "test_count.RDS")
# test <- readRDS("test_table.RDS")
# count <- readRDS("test_count.RDS")
# clust <- cluster(T2, 0.3, "furthest", FALSE, count)

# cluster(phylip=/Users/grejoh/Documents/OptiClusterPackage/Opticluster/sub_set_amazon.txt, method=nearest, cutoff=0.2)
# WritePhylipFile(sparse_matrix@i,sparse_matrix@j, sparse_matrix@x, 0.3, "/Users/grejoh/Documents/mothur/mothur/Clustur/updated_phylip.txt")
# cluster(phylip=/Users/grejoh/Documents/mothur/mothur/Clustur/updated_phylip.txt, method=nearest, cutoff=0.2)