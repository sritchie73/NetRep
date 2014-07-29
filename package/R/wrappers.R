

netStats <- function(discAdj, discIndices, testAdj, testIndices) {
  res <- NetStats(discAdj@address, discIndices, testAdj@address, testIndices)
  lapply(res, as.vector)
}

netProps <- function(adj, subsetIndices) {
  res <- NetProps(adj@address, subsetIndices)
  lapply(res, as.vector)
}

dataProps <- function(discDat, scaledDat, discIndices) {
  res <- DataProps(discDat@address, scaledDat@address, discIndices)
  lapply(res, as.vector)
}
  
