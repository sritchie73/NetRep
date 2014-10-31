#' Wrapper functions for Rcpp functions
#' 
#' These functions abstract away the calls to the Rcpp implementations of 
#' \code{\link{NetStats}}, \code{\link{AdjProps}}, and \code{\link{dataProps}}.
#' There are two reasons why this is convenient:
#' \enumerate{
#'  \item{
#'    \code{\link[bigmemory]{big.matrix}} objects can now be passed directly
#'    without calling their \code{@@address} slot.
#'  }
#'  \item{
#'    The returned armadillo types are all matrices, even though all the results
#'    returned are in vector or scalar format, so we can transform the results 
#'    in R, to avoid risk of memory leaks from incorrect Armadillo to Rcpp
#'    casting on my own part.
#'  }
#' }
#' 
#' @name wrappers
#' 
#' @param discAdj,testAdj,adj \code{\link[bigmemory]{big.matrix}} objects
#'   containing the adjacency matrices for the \emph{discovery} and \emph{test}
#'   networks respectively.
#' @param discIndices,testIndices,subsetIndices indices for the network subset
#'  of interest in the adjacency matrices for the \emph{discovery} and 
#'  \emph{test} and networks respectively.
#' @param datIndices indices for the network subset of interest in the 
#'  underlying data provided for \code{\link{DataProps}}.
#' @param dat underlying data the adjacency matrix was constructed from
#' @param scaledDat a row scaled \code{big.matrix} version of \code{dat}.
#'  See \code{\link{scaleBigMatrix}}.
#' @param lowmem logical; If \code{TRUE}, \code{\link[bigmemory]{big.matrix}} 
#'  objects are attached, and freed at the end of each wrapper function. If
#'  \code{FALSE}, it is assumed the corresponding arguments are the already
#'  attached \code{\link[bigmemory]{big.matrix}} objects.
#'
NULL

#' @rdname wrappers
#' @return \code{netStats:} a list of statistics calculated between adjacency
#'  matrices for a \emph{discovery} and \emph{test} network.
netStats <- function(discAdj, discIndices, testAdj, testIndices, lowmem=FALSE) {
  if (lowmem) {
    discAdj <- attach.big.matrix(discAdj)
    testAdj <- attach.big.matrix(testAdj)
    on.exit({
      rm(discAdj, testAdj)
      gc()
    })
    poke(discAdj, testAdj)
  }
  res <- NetStats(discAdj@address, discIndices, testAdj@address, testIndices)
  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return \code{netProps:} a list of network properties calculated from an
#'  adjacency matrix. These properties can either be scalers (summarising the
#'  whole network subset), or vectors (characterising some property for each
#'  node in the network subset).
adjProps <- function(adj, subsetIndices, lowmem=FALSE) {
  if (lowmem) {
    adj <- attach.big.matrix(adj)
    on.exit({
      rm(adj)
      gc()
    })
    poke(adj)
  }
  res <- AdjProps(adj@address, subsetIndices)
  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return a list of properties quantifying the relationship between a network 
#'   subset and the underlying data the adjacency matrix was calculated from.
#'   These properties can either be scalers (summarising the whole network
#'   subset), or vectors (characterising some property for each node in the
#'   network subset).
dataProps <- function(scaledDat, datIndices, lowmem=FALSE) {
  if (lowmem) {
    scaledDat <- attach.big.matrix(scaledDat)
    on.exit({
      rm(scaledDat)
      gc()
    })
    poke(scaledDat)
  }
  res <- DataProps(scaledDat@address, datIndices)
  lapply(res, as.vector)
}
  
