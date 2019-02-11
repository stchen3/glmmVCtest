#' Calculate the trace of a matrix.
#'
#' @param mat A matrix.
#' @examples
#' tr(diag(5))
tr<-function(mat){
  sum(diag(mat))
}
