#' Calculate network statistics from pre-calculated network properties
#' 
#' @details
#' All network statistics can be fundamentally broken into the calculation of
#' some properties in the discovery network, and some calculation of properties
#' in the test network. It makes sense to split these calculations, so that we
#' only have to calculate the properties in the discovery network/data once.
#' 
#' However, for some statistics it does not make sense to do this due to memory
#' overhead (e.g. \emph{corAdj} would require storing a second copy of most of
#' the discovery network in vector form). For these statistics, see 
#' \code{\link{calcSharedTestStats}}.
#' 
#' @section Background:
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
calcSplitTestStats <- function(discProps, testProps) {
  stopifnot(is.list(discProps) & is.list(testProps))
  stopifnot(length(discProps) == length(testProps))
  
  stats <- c(
    mean.adj = testProps[["mean.adj"]],
    cor.kIM = cor(discProps[["kIM"]], testProps[["kIM"]])
  )
  if ("propVarExpl" %in% names(testProps)) { # Detect if data has been provided
    stats <- c(stats,
      propVarExpl = testProps[["propVarExpl"]],
      mean.kME = mean(sign(discProps[["kME"]]) * testProps[["kME"]]),
      cor.kME = cor(discProps[["kME"]], testProps[["kME"]])
    )
  }
  stats
}

#' Network subset properties
#' 
#' @details
#' All network statistics can be fundamentally broken into the calculation of
#' some properties in the discovery network, and some calculation of properties
#' in the test network. It makes sense to split these calculations, so that we
#' only have to calculate the properties in the discovery network/data once.
#' 
#' However, for some statistics it does not make sense to do this due to memory
#' overhead (e.g. \emph{corAdj} would require storing a second copy of most of
#' the discovery network in vector form). For these statistics, see 
#' \code{\link{calcSharedTestStats}}.
#'
#' @section Background:
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
#' @param adj Adjacency matrix for the network.
#' @param subsetInd Indices of the network subset.
#' @param scaled (Optional) a row scaled \code{big.matrix} of \code{dat}.
#' 
#' @return
#'  A list of topological properties for the given network subset 
#' @seealso \code{\link[=calcSplitTestStats]{Between-network statistics}}
subsetProps <- function(adj, subsetInd, scaled=NULL) {
  props <- adjProps(adj, subsetInd)
  if (!is.null(scaled)) {
    props <- c(props, dataProps(scaled, subsetInd))
  }
  props
}

#' Calculate the cor.cor and mean.cor
#' 
#' For some statistics it does not make sense to calculate the necessary 
#' components in advance due to large memory overhead, or logic that doesn't
#' separate nicely. This function deals with those statistics.
#'  
#' @references 
#' Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my network
#' module preserved and reproducible?} PLoS Comput. Biol. \strong{7}, e1001057
#' (2011).
#' 
#' @seealso \code{\link[=subsetProps]{Network subset topology}} 
#'   \code{\link{netRep}}
#'   
#' @param discCor,testCor \code{\link[bigmemory]{big.matrix}} objects for the 
#'  correlation matrices in the \emph{discovery} and \emph{test} networks 
#'  respectively.
#' @param discIndices,testIndices indices of the network subset in the 
#'  \emph{discovery} and \emph{test} networks respectively.
#' @return A vector of test statistics.
calcSharedTestStats <- function(discCor, discIndices, testCor, testIndices) {
  unlist(netStats(discCor, discIndices, testCor, testIndices))
}
