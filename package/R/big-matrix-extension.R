#' big.matrix diagonals
#' 
#' Extract or replace the diagonal of a \code{\link[bigmemory]{big.matrix}}
#' 
#' @details 
#'  These methods are not supplied by the \code{\link[bigmemory]{bigmemory}}
#'  package. In this package we provide limited support for replacement and 
#'  extraction of the diagonal of \code{big.matrix} objects,
#'  provided they are \code{numeric}. 
#'  
#' @section Additional Notes:
#'  Unlike \code{\link[base]{diag}} for \code{\link[base]{matrix}} objects, you
#'  cannot use \code{diag} to construct a new \code{big.matrix}.
#' 
#' @name big.matrix-extensions
NULL

#' @export
#' @rdname big.matrix-extensions
setMethod("diag<-", signature("big.matrix"), function(x, value) {
  if (typeof(x) != "double") stop("typeof big.matrix must be double")
  SetDiag(x@address, value)
  x
})

#' @export
#' @rdname big.matrix-extensions
setMethod("diag", signature("big.matrix"), function(x) {
  if (typeof(x) != "double") stop("typeof big.matrix must be double")
  GetDiag(x@address)
})
