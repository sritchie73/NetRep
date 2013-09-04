# Added functionality for the S4 class `big.matrix`

library(foreach)

#' @aliases diag<-,big.matrix-method
#' @docType methods
#' @exportMethod diag<-
setMethod("diag<-", signature("big.matrix"), function(x, value) {
  SetDiag(x@address, value)
})

#' @aliases diag,big.matrix-method
#' @docType methods
#' @exportMethod diag
setMethod("diag", signature("big.matrix"), function(x) {
  GetDiag(x@address)
})
