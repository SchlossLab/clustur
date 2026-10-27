#' Read distance matrices
#'
#' Read in distances from a file that is formatted with three columns for the
#' row, column, and distance of a sparse, square matrix or in a phylip-formatted
#' distance matrix.
#'
#' @export
#' @param distance_file Either a phylip or column distance file, or a sparse
#' matrix. The function will detect the format for you.
#' @param count_table A table of names and the given abundance per group. Can
#' be in mothur's sparse or full format. The function will detect the format for
#' you.
#' @param cutoff The value you wish to use as a cutoff when clustering.
#' @param is_similarity_matrix are you using a similarity matrix (default) or
#' distance matrix?
#' @return A distance `externalptr` object that contains all your
#' distance information. Can be accessed using `get_distance_df()`
#'
#' @examples
#'
#'  i_values <- as.integer(1:100)
#'  j_values <- as.integer(sample(1:100, 100, TRUE))
#'  x_values <- as.numeric(runif(100, 0, 1))
#'  s_matrix <- create_sparse_matrix(i_values, j_values, x_values)
#'
#'  sparse_count <- data.frame(
#'                  Representative_Sequence = 1:100,
#'                  total = rep(1, times = 100))
#'
#'  column_path <- example_path("amazon_column.dist")
#'  phylip_path <- example_path("amazon_phylip.dist")
#'  count_table <- read_count(example_path("amazon.full.count_table"))
#'
#'  data_column <- read_dist(column_path, count_table, 0.03)
#'  data_phylip <- read_dist(phylip_path, count_table, 0.03)
#'  data_sparse <- read_dist(s_matrix, sparse_count, 0.03)

read_dist <- function(distance_file, count_table,
                      cutoff, is_similarity_matrix = FALSE) {
  count_table <- validate_count_table(count_table)
  if (inherits(distance_file, "character")) {
    if (!file.exists(distance_file)) {
      stop("Invalid file path: please enter a new file path.")
    }
    result <- ProcessDistanceFiles(distance_file,
                                   count_table, cutoff, is_similarity_matrix)
    class(result) <- "distance_object"
    return(result)
  }

  # Its a sparse matrix not a path
  # filter out the sparse matrix
  if (!inherits(distance_file, "dgTMatrix")) {
    stop("If you are not using a file, ensure the sparse matrix is created from
     the `create_sparse_matrix()` function")
  }
  result <- ProcessSparseMatrix(distance_file@i, distance_file@j,
                                distance_file@x,
                                count_table, cutoff, is_similarity_matrix)
  class(result) <- "distance_object"
  result
}


#' Cluster entities together
#'
#' Clusters entities represented in a distance matrix and count table using
#' one of several algorithms and outputs information about the composition and
#' abundance of each cluster
#'
#' @export
#' @param distance_object The distance object that
#'  was created using the `read_dist()` function.
#' @param cutoff The cutoff you want to cluster towards.
#' @param method The method of clustering to be performed: opticlust (default),
#' furthest, nearest, average, or weighted.
#' @param feature_column_name_to Set the name of the column in the cluster
#' dataframe that contains the sequence names.
#' @param bin_column_name_to Set the name of the column in the
#' cluster dataframe that contains the name of the group of sequence names.
#' @param random_seed the random seed to use, (default = 123).
#' @param delta Allows you to set the stable value for the metric in the
#' opticluster method. It defaults to 0.0001, to reach full convergence
#' set this value to 0.
#' @param max_iterations the max amount of iterations the opticluster
#' method will run. The default is 100.
#' @param initialize this parameter allows you to set the default randomization of
#' for the opticluster method. The two methods are: "singleton" and "oneotu." Singleton
#' will randomly assign each sequence to its own OTU while "oneotu" assigns all sequence
#' into one OTU. Defaults to singleton.
#' @param strollur_object Using the `strollur` r package, you can use your
#' created strollur object to save your clustering results to. This allows you
#' to seamless modify your strollur object and analysis results as needed. 
#' @return A list of `data.frames` that contain abundance, and clustering
#' results. If you used `method = opticlust`, it will also return clustering
#' performance metrics.
#'
#' @examples
#'
#'  cutoff <- 0.03
#'  count_table <- read_count(example_path("amazon.full.count_table"))
#'  distance_data <- read_dist(example_path("amazon_column.dist"),
#'                             count_table, cutoff)
#'
#'  cluster_results <- cluster(distance_data,
#'                             cutoff, method = "opticlust",
#'                             feature_column_name_to = "sequence",
#'                             bin_column_name_to = "omu")
#'  cluster_results <- cluster(distance_data,
#'                             cutoff, method = "furthest")
#'  cluster_results <- cluster(distance_data,
#'                             cutoff, method = "nearest")
#'  cluster_results <- cluster(distance_data,
#'                             cutoff, method = "average")
#'  cluster_results <- cluster(distance_data,
#'                             cutoff, method = "weighted")
#'
#'
cluster <- function(distance_object, cutoff, method = "opticlust",
                    feature_column_name_to = "feature",
                    bin_column_name_to = "bin", random_seed = 123,
                    delta = 0.0001, max_iterations = 100, 
                    initialize = "singleton", strollur_object = NULL) {
  if (!inherits(distance_object, "distance_object")) {
    stop("`distance_object` must be generated using the `read_dist` function")
  }
  if (!(method %in% c("opticlust",
                      "furthest", "nearest", "average", "weighted", "optifit"))) {
    stop("`method` parameter can only be opticlust, furthest, nearest,
         average, or weighted.")
  }
  df <- data.frame()
 if (method != "opticlust") {
    df <- Cluster(distance_object, method,
                  feature_column_name_to, bin_column_name_to, cutoff,
                  random_seed)
  } else {
    df <- OptiClust(distance_object, feature_column_name_to,
                      bin_column_name_to, cutoff, random_seed,
                      delta, max_iterations, initialize)
    df$iteration_metrics <- df$iteration_metrics[, c("iter", "time",
                                                     "label", "num_otus",
                                                     "cutoff", "tp", "tn",
                                                     "fp", "fn", "sensitivity",
                                                     "specificity", "ppv",
                                                     "npv", "fdr", "accuracy",
                                                     "mcc", "f1score")]

  }
  class(df) <- "mothur_cluster"
  if (!is.null(strollur_object)) {
    result <- add_results_to_strollur_object(df)
    strollur::xdev_assign_bin_representative_sequences(strollur_object,
                                                       result)                                               
  }
  df
}

add_results_to_strollur_object <- function(cluster_result, strollur_object) {
  results <- split_clusters_to_list(cluster_result)
  list_names <- names(results)
  result <- list()
  for (i in seq_along(results)) {
    result <- append(result, lapply(results[[i]], function(x) {
              c(x, list_names[[i]])
            }))
  }
  result <- Reduce(rbind, result)
  result <- data.frame(bin_name = result[, 2], sequence_name = result[, 1])
  result
}

#' Validate Count Table
#'
#' Determines whether user supplied count table is valid
#'
#' @export
#' @param count_table_df The count table `data.frame` object.
#' @return A validated count table `data.frame`
#' @description
#' If the count table is already valid nothing will change,
#' otherwise it will add a new group to the count table file.
#'
#' @examples
#'  count_table <- read.delim(example_path("amazon.full.count_table"))
#'  count_table_valid <- validate_count_table(count_table)
#'
validate_count_table <- function(count_table_df) {
  if (ncol(count_table_df) > 2) {
    return(count_table_df)
  }
  totals <- count_table_df$total
  count_table_df <- cbind(count_table_df, totals)
  names(count_table_df)[3] <- "no_group"
  count_table_df[[1]] <- as.character(count_table_df[[1]])
  count_table_df
}

#' Example Path
#'
#' This function was created as a helper function to generate file paths to our
#' internal data. You should use this function if you
#' want to follow along with the example, or interact with the data
#' @param file The file name of the data; leave as NULL (default) to get full
#' list of example files
#' @examples
#' example_path("amazon_phylip.dist")
#' example_path()
#' @return the path to the file as a `character` or a vector of `character`
#' giving example filenames if `fill = NULL`.
#' @export

example_path <- function(file = NULL) {
  path <- ""
  if (is.null(file)) {
    path <- dir(system.file("extdata", package = "clustur"))
  } else {
    path <- system.file("extdata", file, package = "clustur", mustWork = TRUE)
  }
  path
}

#' Create Sparse Matrix
#'
#' Given a list of i indexes, j indexes, and distances values,
#' we can create a sparse distance matrix for you. Each vector
#' must have the same size.
#' @param i_index A list of i indexes, must be numeric
#' @param j_index A list of j indexes, must be numeric
#' @param distances A list of the distance at the i and j index
#' @examples
#'
#'  i_values <- as.integer(1:100)
#'  j_values <- as.integer(sample(1:100, 100, TRUE))
#'  x_values <- as.numeric(runif(100, 0, 1))
#'  s_matrix <- create_sparse_matrix(i_values, j_values, x_values)
#' @return a `dgTMatrix` from the `Matrix` library.
#' @export
create_sparse_matrix <- function(i_index, j_index, distances) {
  size <- max(i_index, j_index)
  spMatrix(size, size, i_index, j_index, distances)
}

#' @title Optifit
#' 
#' @description Optifit is the 'cluster.fit' command in mothur. It allows you to fit sequences to OTUs or existing OTUS. Just like
#' mothur, there are two ways of running this command: open, or closed. Closed fits reads into existing OTUs, 
#' and scraps any reads that cannot be fitted(default it closed). Open fits reads into existing OTUs while unfitted
#' OTUs are clustered into new OTUs.
#' @param fit_dist The distance object that will be used to fit to the new reference. This must be
#' created using the `read_dist()` function.
#' @param cutoff The cutoff you want to cluster towards.
#' @param feature_column_name Set the name of the column in the cluster
#' dataframe that contains the sequence names. This is a `character` value.
#' @param bin_column_name Set the name of the column in the
#' cluster dataframe that contains the name of the group of sequence names. This is a `character` value.
#' @param accnos a `vector` of strings that represent the reference sequences by name. Use this when you know the
#' group of sequences you want your fit data to fit against.
#' @param ref_fasta the file name of your reference's fasta data. This is a `character` value.
#' @param ref_count the file name of your reference's count table. This is a `character` value.
#' @param ref_list the file name of your reference's list file that you want the fit data to fit against.
#' This is a `character` value.
#' @param fit_fasta the file name of your fitted data's fasta file. This is a `character` value.
#' @param fit_percent allows you to set the percentage of reads that need to be fitted. This is a `numerical` value.
#' @param closed deteremines whether you want the function to run using the closed or open method. The closed method
#' will return reads that fit into existing otus and scrap those that did not. The open method will keep reads
#' that fit into existing otus and create new OTUs for those that did not. This is a `logical` value.
#' @param print_ref controls whether all the reference sequences are printed in the final result.
#' This is a `character` value.
#' @param random_seed a `integer` value that represents the current random seed.
#' @param max_iterations the max number of iterations that cluster computes. This is a `integer` value.
#' @param denovo_iters the max number of iterations the denovo method computes.  This is a `integer` value.
#' @param optimal_denovo_result determines whether or not you want the denovo method to return the cluster
#' with the highest mcc value or all of the randomly computed clusters. This is a `logical` value.
#' @param number_of_threads a `integer` value that represents the amount of allotted threads that
#' can be used.
#' @examples
#'  amazon_count <- read_count(example_path("amazon.full.count_table"))
#'  amazon_dist <- read_dist(example_path("amazon_phylip.dist"), amazon_count, 0.03)
#'
#'  esophagus_fasta <- example_path("esophagus.fasta")
#'  esophagus_count <- read_count(example_path("esophagus.count_table"))
#'  amazon_fasta <- example_path("amazon.align")
#'  esophagus_list <- example_path("esophagus.list")
#'
#'
#'  denovo <- optifit(amazon_dist, 0.03)
#'  denovo_open <- optifit(amazon_dist, 0.03, close = FALSE)
#'  accnos <- optifit(amazon_dist, 0.03, accnos = amazon_count$Representative_Sequence[1:40])
#'  accnos_open <- optifit(amazon_dist, 0.03, accnos = amazon_count$Representative_Sequence[1:40],
#'                         closed=FALSE)
#'
#'  ref_cluster <- optifit(amazon_dist, 0.03, "feature", "bin", accnos = NULL,
#'                         ref_fasta = esophagus_fasta, ref_count = esophagus_count,
#'                         ref_list = esophagus_list, fit_fasta = amazon_fasta,
#'                         number_of_threads = parallel::detectCores(), print_ref = TRUE)
#' @return A `list` of a list of `data.frames` that contain abundance, and clustering
#' results. If you used `method = opticlust`, it will also return clustering
#' performance metrics.
#' @export
optifit <- function(fit_dist, cutoff, feature_column_name = "feature", bin_column_name = "bin",  accnos = NULL, ref_fasta = NULL,
                    ref_count = NULL, ref_list = NULL, fit_fasta = NULL, fit_percent = 50, closed = TRUE, print_ref = FALSE,
                    random_seed = 123, max_iterations = 100, denovo_iters = 10, optimal_denovo_result = TRUE, 
                    number_of_threads = 1) {

  result <- c()
  ref_parameters_check <- all( as.logical(lapply(list(ref_fasta, ref_list, fit_fasta, ref_count), is.null)))
  if (ref_parameters_check && is.null(accnos)) {
     message(paste0("accnos, ref_fasta, ref_count, ref_list, and fit_fasta parameters are",
                   " all NULL. Running Optifit using denovo."))
    if(fit_percent < 0.01) {
      fit_percent = 0.01
    }
    else if(fit_percent > 100) {
      fit_percent = 100
    }
    # run denovo
    result <- OptiFit(fit_dist, feature_column_name, bin_column_name, cutoff, fit_percent, closed, print_ref, 
                      selfReference = FALSE, max_iterations, denovo_iters, random_seed)
  } else if (!is.null(accnos)) {
    message("Accnos file is not null, running optifit using accnos.")
    result <- OptiFit2(fit_dist, feature_column_name, bin_column_name, accnos, cutoff, closed, print_ref,
                       selfReference = FALSE, max_iterations, random_seed)
  } else if (!ref_parameters_check) {
    message("ref_fasta, ref_list, ref_count, and fit_fasta are not NULL, running reference based optifit.")
    # all ref list names have to be present in the fasta error catch that later
    ref_fasta <- strollur::read_fasta(ref_fasta)
    fit_fasta <- strollur::read_fasta(fit_fasta)
    ref_list <- strollur::read_mothur_list(ref_list)
    copied_data <- CopyObject(fit_dist)
    class(copied_data) <- "distance_object"
    AddDataToDistanceData(copied_data, ref_list, ref_count, ref_fasta, fit_fasta, cutoff, number_of_threads)
    result <- OptiFit3(copied_data, ref_list, ref_fasta$sequence_name, feature_column_name, bin_column_name, cutoff, closed,
                       print_ref, selfReference = TRUE, max_iterations, random_seed)
  }
  else {
    stop(paste0("You either have forget to supply proper parameters are you did not supply all,",
                "the parameters required to run the reference clustering method, please check again."))
  }
  class(result) <- "mothur_cluster"
  result
}

#' @export
#' @title Cluster Split
#' @description
#' Cluster Split splits data into groups based on taxomic group. It then clusters them using 
#' the following methods: opticluster, average neighbor, nearest neighbor, furthest
#' neighbor, and weighted neighbor. You can read more about
#' [cluster split](https://mothur.org/wiki/cluster.split/) on the mothur wiki.
#' @param fasta_file a `string` value that represents the path to your fasta file.
#' @param taxonomy_file a `string` value that represents the path to your taxonomy file.
#' @param count_table the count table that represents the data, must be read using the `read_count()`
#' function.
#' @param cluster_method the clustering method you wish to use. Your clustering options are: 
#' opticluster, average neighbor, nearest neighbor, furthest
#' neighbor, and weighted neighbor.
#' @param cutoff a `numerical` value that represents the distance cutoff.
#' @param taxonomy_cutoff a `integer` value that represents the taxonomic level cutoff. 
#' @param feature_column_name Set the name of the column in the cluster
#' dataframe that contains the sequence names. This is a `character` value.
#' @param bin_column_name Set the name of the column in the
#' cluster dataframe that contains the name of the group of sequence names. This is a `character` value.
#' @param random_seed a `integer` value that represents the current random seed.
#' @param number_of_threads a `integer` value that represents the amount of allotted threads that
#' can be used.
#' @examples
#' 
#'  cutoff <- 0.2
#'  count_table <- read_count(example_path("amazon.full.count_table"))
#'  fasta <- example_path("amazon.align")
#'  taxonomy <- example_path("amazon.taxonomy")
#'  df <- cluster_split(fasta, taxonomy, count_table, "opti", cutoff, 3, number_of_threads = 1)
#' 
#' @return A list of `data.frames` that contain abundance, and clustering
#' results.
cluster_split <- function(fasta_file, taxonomy_file, count_table, cluster_method, cutoff,
                       taxonomy_cutoff, feature_column_name = "feature", bin_column_name = "bin",  
                       random_seed = 123, number_of_threads = 1) {
  

  fasta <- strollur::read_fasta(fasta_file)
  taxonomy <- strollur::read_mothur_taxonomy(taxonomy_file)

  OptiSplit(fasta, taxonomy, count_table, cluster_method, feature_column_name,
            bin_column_name, cutoff, taxonomy_cutoff, random_seed, number_of_threads)
}
