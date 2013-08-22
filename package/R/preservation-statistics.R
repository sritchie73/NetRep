# This file holds all the statistics to be calculated between networks
# 
# Note: come back and formally attribute the WGCNA package and authors.
# 

#' Calculate the meanAdj statistic of a module in a test network
#' 
#' Calculates the mean edge weight of the module in the test network.
#' A high preservation score indicates preservation of the module's density.
#' 
#' @note Make sure to mark with NA edges that should not be included in the 
#' calculation. This is done to minimise computation.
#' For example, gene coexpression adjacency matrices from WGCNA should have 
#' their diagonal replaced with NA. For sparse networks the user will need to 
#' determine where NA or 0 is most appropriate.
#' 
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'    Steve. Is My Network Module Preserved and Reproducible?. PLoS 
#'    Computational Biology, 2011.
#' 
#' @param testNet the adjacency matrix for the test network.
#' @param testIndices the indices on the nodes in the module.
#' @return The mean edge weight of the module in the test network.
meanAdj <- function(testNet, testIndices) {
  mean(testNet[testIndices, testIndices], na.rm=TRUE)
}
