# Added functionality for the S4 class `big.matrix`

setGeneric("diag<-")
setMethod("diag<-", signature("big.matrix"), function(x, value) {
  foreach(i = 1:ncol(x)) %do%{
    x[i,i] <- value
  }
  x
})

setGeneric("diag")
setMethod("diag", signature("big.matrix"), function(x) {
  foreach(i = 1:ncol(x), .combine=c) %do%{
    x[i,i] <- value
  }
})