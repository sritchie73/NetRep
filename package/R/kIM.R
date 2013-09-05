#' Calculates the connectivity (kIM) of a module.
#' 
#' Calculates the inter or intra modular connectivity of all nodes in a given
#' module. This is used for calculating several preservation statistics.
#' 
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'    Steve. Is My Network Module Preserved and Reproducible?. PLoS 
#'    Computational Biology, 2011.
#' 
#' @param adjacency Adjacency matrix to calculate meanAdj of module in.
#' @param moduleIndices Where the nodes in the module are located in the 
#'  adjacency network.
#' @param allNodes Logical; if FALSE, the connectivity is calculated based only
#'  on nodes in the module, otherwise a module's connectivity to all nodes is 
#'  calculated.
#' @return The connectivity for each node in the module.
kIM <- function(adjacency, moduleIndices, allNodes) {
  KIM(adjacency@address, moduleIndices, allNodes)
}

kIMR <- function(adjacency, moduleIndices, allNodes) {
  # For testing equivalence in R.
  if (allNodes) {
    subset <- adjacency[,moduleIndices]
  } else {
    subset <- adjacency[moduleIndices, moduleIndices]
  }
  colSums(subset, na.rm=TRUE)
}
