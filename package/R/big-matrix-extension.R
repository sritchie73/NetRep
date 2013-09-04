# Added functionality for the S4 class `big.matrix`

library(foreach)

setMethod("diag<-", signature("big.matrix"), function(x, value) {
  SetDiag(x@address, value)
})

setMethod("diag", signature("big.matrix"), function(x) {
  GetDiag(x@address)
})
