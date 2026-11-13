#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Matrix spMatrix
#' @importFrom methods el
#' @importFrom Rcpp evalCpp
#' @importFrom reshape2 dcast
#' @importFrom utils read.delim
#' @useDynLib clustur
## usethis namespace: end
NULL

#' @importFrom strollur read_mothur_taxonomy
#' @export
strollur::read_mothur_taxonomy

#' @importFrom strollur read_fasta
#' @export
strollur::read_fasta

#' @importFrom strollur read_mothur_list
#' @export
strollur::read_mothur_list
