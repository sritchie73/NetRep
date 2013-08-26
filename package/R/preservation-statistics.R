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
#' @param adjacency Adjacency matrix to calculate meanAdj of module in.
#' @param moduleIndicies Where the nodes in the module are located in the 
#'  adjacency network.
#' @param includeDiagonal Logical, if FALSE, diagonal values (edges between 
#'  nodes and themselves) are not counted towards the mean edge weight.
#' @return The mean edge weight of the module in the test network.
meanAdj <- function(adjacency, moduleIndices, includeDiagonal) {
  MeanAdj(adjacency@address, moduleIndices, includeDiagonal)
}
