#' big.matrix diagonals
#' 
#' Extract or replace the diagonal of a \code{\link[bigmemory]{big.matrix}}
#'  
#' @section Additional Notes:
#'  Unlike \code{\link[base]{diag}} for \code{\link[base]{matrix}} objects, you
#'  cannot use \code{diag} to construct a new \code{big.matrix}.
#' 
#' @name big.matrix-extensions
#' @import methods
NULL

#' @export
#' @rdname big.matrix-extensions
setMethod("diag<-", signature("big.matrix"), function(x, value) {
  SetDiag(x@address, value)
  x
})

#' @export
#' @rdname big.matrix-extensions
setMethod("diag", signature("big.matrix"), function(x) {
  GetDiag(x@address)
})
