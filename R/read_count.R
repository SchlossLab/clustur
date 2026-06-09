#' Read count table
#'
#' This function will read and return your count table. It can take in
#' sparse and full count tables.
#'
#' @param count_table The file path of your count table or strollur object.
#' @examples
#' count_table <- read_count(example_path("amazon.full.count_table"))
#' @return a count table `data.frame`.
#' @export
read_count <- function(count_table) {
  UseMethod("read_count", count_table)
}

#' @export
#' @rdname read_count
read_count.strollur <- function(count_table) {
  abundances_table <- strollur::abundance(count_table, by_sample = T)
  shared_table <- reshape2::dcast(abundances_table, sequence_name ~ sample, value.var = "abundance")
  shared_table[is.na(shared_table)] <- 0
  row_names <- strollur::names(count_table, type = "sample")
  shared_table$total <- rowSums(shared_table[, row_names])
  cbind(data.frame(Representative_Sequence = shared_table$sequence_name,
                  total = shared_table$total),
        shared_table[, row_names])
}

#' @export
#' @rdname read_count
read_count.character <- function(count_table) {
  test_read <- read.delim(count_table, sep = "\t", header = FALSE)
  if (grepl("#", test_read[1, 1], fixed = TRUE)) {
    count_table_sparse <- read.delim(count_table, sep = "\t", skip = 2)
    count_table_sparse <- lapply(count_table_sparse, as.character)
    ct <- CreateDataFrameFromSparseCountTable(count_table_sparse)
    return(validate_count_table(ct))
  }
  validate_count_table(read.delim(count_table, sep = "\t"))
}