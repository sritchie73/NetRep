# This file holds all the statistics to be calculated between networks
# 
# Note: come back and formally attribute the WGCNA package and authors.
# 

#' Calculate the meanAdj statistic of a module in a test network
#' 
#' Calculates the mean edge weight of the module in the test network.
#' A high preservation score indicates preservation of the module's density.
#' 
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'    Steve. Is My Network Module Preserved and Reproducible?. PLoS 
#'    Computational Biology, 2011.
#' 
#' @param testNet the adjacency matrix for the test network.
#' @param testIndices the indices on the nodes in the module.
#' @param diagonal How to handle diagonals in the adjacency matrix:
#'   \itemize{
#'     \item{NA }{An edge from a node to itself is allowed and counted in the
#'                mean.}
#'     \item{numeric }{The value the all diagonal entries of the adjacency
#'                     matrix have, to be removed in the calculation of the 
#'                     mean.}
#'   }
#' @return The mean edge weight of the module in the test network.
meanAdj <- function(testNet, testIndices, diagonal=NA) {
  nNodes <- length(testIndices)
  if(is.na(diagonal)) {
    mean <- sum(testNet[testIndices, testIndices])/{nNodes * nNodes}
  } else {
    mean <- sum(testNet[testIndices, testIndices]) - nNodes*diagonal
    mean <- mean / {nNodes*nNodes - nNodes}
  }
  return(mean)
}
