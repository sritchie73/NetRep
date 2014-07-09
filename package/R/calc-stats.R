#' Calculate between-dataset network statistics.
#' 
#' For a given network subset, as defined by the nodes in the \code{nodeDict},
#' calculate 
#'
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'    Steve. Is My Network Module Preserved and Reproducible?. PLoS 
#'    Computational Biology, 2011.
#' 
#' @param discDat Underlying data for the reference network.
#' @param discAdj Adjacency matrix for the discovery network.
#' @param testDat Underlying data for the test network.
#' @param testAdj Adjacency matrix for the test network.
#' @param nodeDict A list with four elements, which map the nodes of the network
#'   subset of interest to their indices in the \code{discDat}, \code{discAdj},
#'   \code{testDat}, and \code{testAdj} respectively.
#' @return A vector of preservation scores for the network subset.
#'   
calcReplStats <- function(discDat, discAdj, testDat, testAdj,
                          nodeDict) {
  stopifnot(class(discDat) %in% c("big.matrix", "NULL"))
  stopifnot(class(discAdj) %in% c("big.matrix", "NULL"))
  stopifnot(class(testDat) %in% c("big.matrix", "NULL"))
  stopifnot(class(testAdj) %in% c("big.matrix", "NULL"))
  
  # Handle on the fly calculations of the adjacency matrix
  if (is.null(testAdj)) {
    stop("not yet implemented")
  }
  if (is.null(discAdj)) {
    stop("not yet implemented")
  }
  
  # Calculate the statistics for both networks that use only the Adjacency
  # Matrices
  preservation <- c(
      meanAdj(testAdj, nodeDict[["testNet"]]),
      cor(
          kIM(discAdj, nodeDict[["discNet"]], FALSE), 
          kIM(testAdj, nodeDict[["testNet"]], FALSE), 
          method="spearman"
        ),
      cor(
          kIM(discAdj, nodeDict[["discNet"]], TRUE), 
          kIM(testAdj, nodeDict[["testNet"]], TRUE), 
          method="spearman"
        )
    )
  names(preservation) <- c("meanAdj", "corKIM", "corKALL")
  
  # Calculate the statistics that also require the underlying data the netwrok w
  # was constructed from.
  if (!is.null(discDat) & !is.null(testDat)){ 
    datPres <- c(
    )
    names(datPres) <- c()
    preservation <- c(preservation, datPres)
  } 
  
  return(preservation)
}