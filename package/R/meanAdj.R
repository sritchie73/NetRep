#' Calculate the meanAdj statistic of a module in a test network
#' 
#' Calculates the mean edge weight of the module in the test network.
#' A high preservation score indicates preservation of the module's density.
#' 
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'    Steve. Is My Network Module Preserved and Reproducible?. PLoS 
#'    Computational Biology, 2011.
#' 
#' @param adjacency Adjacency matrix to calculate meanAdj of module in.
#' @param moduleIndices Where the nodes in the module are located in the 
#'  adjacency network.
#' @return The mean edge weight of the module in the test network.
meanAdj <- function(adjacency, moduleIndices) {
  # Error check inputs before passing to C++
  stopifnot(class(adjacency) == "big.matrix")
  stopifnot(is.vector(moduleIndices))
  stopifnot(class(moduleIndices) %in% c("numeric", "integer"))
  
  MeanAdj(adjacency@address, moduleIndices)
}

meanAdjR <- function(adjacency, moduleIndices) {
  # R implementation of meanAdj. Used for testing correctness of C++
  # implementation. Does not handle missingness unless diagonal included!
  mean(adjacency[moduleIndices, moduleIndices], na.rm=TRUE)
}