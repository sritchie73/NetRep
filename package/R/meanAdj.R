#' Mean Adjacency of a Network
#' 
#' The mean adjacency statistic measures how \emph{dense} a network is. It is
#' simply the average edge weight of the given network, or network subset.
#' 
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'    Steve. Is My Network Module Preserved and Reproducible?. PLoS 
#'    Computational Biology, 2011.
#' 
#' @param adjacency Adjacency matrix.
#' @param subsetIndices row/column indices in the provided \code{adjacency} that
#'   correspond to the network subset of interest.
#' @param undirected Logical; if \code{TRUE}, only the lower half of the 
#'   symmetric adjacency matrix will be used in the calculations.
#' @return a single numeric value.
#' @export
meanAdj <- function(adjacency, subsetIndices, undirected) {
  stopifnot(class(adjacency) == "big.matrix")
  MeanAdj(adjacency@address, subsetIndices, undirected)
}

#' Mean Adjacency of a Network, R Implementation.
#' 
#' Used to unit test the C++ functionality of \code{\link{meanAdj}}.
#' 
#' @param adjacency Adjacency matrix.
#' @param subsetIndices row/column indices in the provided \code{adjacency} that
#'   correspond to the network subset of interest.
#' @return a single numeric value.
meanAdjR <- function(adjacency, subsetIndices) {
  mean(abs(adjacency[subsetIndices, subsetIndices]), na.rm=TRUE)
}
