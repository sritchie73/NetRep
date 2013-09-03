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
#' @param includeDiagonal Logical, if FALSE, diagonal values (edges between 
#'  nodes and themselves) are not counted towards the mean edge weight.
#' @return The mean edge weight of the module in the test network.
meanAdj <- function(adjacency, moduleIndices, includeDiagonal) {
  # Error check inputs before passing to C++
  stopifnot(class(adjacency) == "big.matrix")
  stopifnot(is.vector(moduleIndices))
  stopifnot(class(moduleIndices) %in% c("numeric", "integer"))
  stopifnot(class(includeDiagonal) == "logical")
  stopifnot(length(includeDiagonal) == 1)
  stopifnot(!is.na(includeDiagonal))
  
  MeanAdj(adjacency@address, moduleIndices, includeDiagonal)
}

meanAdjR <- function(adjacency, moduleIndices, includeDiagonal) {
  # R implementation of meanAdj. Used for testing correctness of C++
  # implementation. Does not handle missingness unless diagonal included!
  if(!includeDiagonal) {
    meanAdj <- (sum(adjacency) - sum(diag(adjacency))) / (n*n - n)
  } else {
    meanAdj <- mean(adjacency[moduleIndices, moduleIndices], na.rm=TRUE)
  }
  return(meanAdj)
}