#' Node Connectivity (weighted degree)
#' 
#' If \code{allNodes} is \code{FALSE}, then the \emph{intramodular connectivity}
#' (kIM) is calculated: this is the sum of the edgeweights for each node in the 
#' network subset, to all other nodes in that network subset. If instead we
#' specify \code{allNodes = TRUE}, then the \emph{whole network connectivity}
#' for those nodes is calculated.
#' 
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'   Steve. Is My Network Module Preserved and Reproducible?. PLoS Computational
#'   Biology, 2011.
#'   
#' @param adjacency Adjacency matrix representation of the network
#' @param subsetIndices row/column indices in the provided \code{adjacency} that
#'   correspond to the network subset of interest.
#' @param allNodes Logical; if \code{FALSE}, the connectivity is calculated
#'   based only on nodes in the module, otherwise a module's connectivity to all
#'   nodes is calculated.
#' @param undirected Logical; if \code{TRUE}, only the lower half of the 
#'   symmetric adjacency matrix will be used in the calculations.
#' @return The connectivity (weighted degree) for each node in the network
#' subset
#' @export
kIM <- function(adjacency, subsetIndices, allNodes, undirected) {
  if (allNodes) {
    # marginally faster than a C++ implementation.
    sums <- colSums(abs(adjacency[,subsetIndices]), na.rm=TRUE)
    if (undirected) {
      sums <- sums/2
    }
    return(sums) 
  } else { 
    return(KIM(adjacency@address, subsetIndices, undirected))
  }
}

#' Node Connectivity (weighted degree), R Implementation
#' 
#' Used to unit test the C++ functionality of \code{\link{kIM}}.
#' 
#' @param adjacency Adjacency matrix representation of the network
#' @param subsetIndices row/column indices in the provided \code{adjacency} that
#'   correspond to the network subset of interest.
#' @param allNodes Logical; if \code{FALSE}, the connectivity is calculated
#'   based only on nodes in the module, otherwise a module's connectivity to all
#'   nodes is calculated.
#' @param undirected Logical; if \code{TRUE}, only the lower half of the 
#'   symmetric adjacency matrix will be used in the calculations.
#' @return The connectivity (weighted degree) for each node in the network
#' subset
kIMR <- function(adjacency, subsetIndices, allNodes, undirected) {
  # For testing equivalence in R.
  if (allNodes) {
    subset <- adjacency[,subsetIndices]
  } else {
    subset <- adjacency[subsetIndices, subsetIndices]
  }
  sums <- colSums(abs(subset), na.rm=TRUE)
  if (undirected) {
    sums <- sums/2
  }
  sums
}

