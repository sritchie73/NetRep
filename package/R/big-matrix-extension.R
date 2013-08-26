# Added functionality for the S4 class `big.matrix`

library(foreach)

setMethod("diag<-", signature("big.matrix"), function(x, value) {
  foreach (i = 1:ncol(x)) %do% {
    x[i,i] <- value
  }
  x
})

setMethod("diag", signature("big.matrix"), function(x) {
  foreach (i = 1:ncol(x), .combine=c) %do% {
    x[i,i]
  }
})
