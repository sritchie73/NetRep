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
#' @return a single numeric value.
#' @export
meanAdj <- function(adjacency, subsetIndices) {
  stopifnot(class(adjacency) == "big.matrix")
  MeanAdj(adjacency@address, subsetIndices)
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
  mean(adjacency[subsetIndices, subsetIndices], na.rm=TRUE)
}
