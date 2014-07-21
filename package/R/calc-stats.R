#' Between-dataset network statistics.
#' 
#' Calculate network subset replication test statistics from their network 
#' properties (as calculated by \code{\link{subsetProps}}).
#' 
#' @details 
#' The returned test statistics indicate the replication/preservation of each 
#' network property for a network subset in the two datasets provided in the 
#' \code{discProps} and \code{testProps} arguments.
#' 
#' Evaluating the significance of these test statistics requires generation of a
#' null distribution for each, by randomly sampling network subsets of the same
#' size in the \emph{test} dataset, calculating their replication, and 
#' evaluating how extreme the observed test statistic is, in comparison to test 
#' statistics drawn from random permutation. This procedure is performed by the 
#' main function of this package, \code{\link{netRep}}.
#' 
#' @references 
#' Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my network
#' module preserved and reproducible?} PLoS Comput. Biol. \strong{7}, e1001057
#' (2011).
#' 
#' @seealso \code{\link[=subsetProps]{Network subset topology}} 
#'   \code{\link{netRep}}
#' 
#' @param discProps,testProps properties of the network subset in the each of 
#'  the \emph{discovery} and \emph{test} datasets, respectively, as calculated 
#'  by \code{\link{subsetProps}}.
#' @return A vector of test statistics.
#'   
#' @import foreach
#' @export
subsetTestStats <- function(discProps, testProps) {
  stopifnot(is.list(discProps) & is.list(testProps))
  stopifnot(length(discProps) == length(testProps))
  
  stats <- c(
    meanAdj = testProps[["meanAdj"]],
    cor.kIM = cor(discProps[["kIM"]], testProps[["kIM"]], method="pearson"),
    cor.kALL = cor(discProps[["kALL"]], testProps[["kALL"]], method="pearson")
  )
}

#' Network subset topology
#' 
#' Calculate the topological characteristics of a network subset.
#'  
#' @references 
#'  Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my network
#'  module preserved and reproducible?} PLoS Comput. Biol. \strong{7}, e1001057
#'  (2011).
#'  
#' @param adj Adjacency matrix for the network.
#' @param adjInd Indices of the network subset in \code{adj}
#' @param dat (Optional) Underlying data for the network.
#' @param datInd (Optional) Indices of the network subset in \code{dat}
#' @param undirected logical; If \code{TRUE}, only the lower half of \code{adj}
#'   is used to calculate some propertiest (where applicable).
#'
#' @return
#'  A list of topological properties for the given network subset 
#' @seealso \code{\link[=subsetTestStats]{Between-network statistics}}
#' @export
subsetProps <- function(adj, adjInd, dat=NULL, datInd=NULL, undirected=FALSE) {
  # Sanity check user input.
  stopifnot(class(adj) %in% c("big.matrix"))
  stopifnot(is.vector(adjInd) & class(adjInd) %in% c("integer", "numeric"))
  stopifnot(class(dat) %in% c("big.matrix", "NULL"))
  if (is.null(dat)) {
    if (length(datInd) > 1) {
      warning("datInd provided, but dat is NULL, ignoring.")
    }
  }
  if (length(datInd) == 0) {
    if(!is.null(dat)) {
      stop("dat provided, but no indices for the network subset provided.",
           " aborting.")
    }
  } else {
    stopifnot(is.vector(datInd) & class(datInd) %in% c("integer", "numeric"))
  }
  
  # TODO: on the fly network construction.
  
  props <- list(
    meanAdj = meanAdj(adj, adjInd, undirected),
    kIM = kIM(adj, adjInd, FALSE),
    kALL = kIM(adj, adjInd, TRUE)
  )
  if (!missing(dat)) {
    # TODO:
  }
  props
}
