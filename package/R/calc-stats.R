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
    meanKIM = testProps[["meanKIM"]],
    meanMAR = testProps[["meanMAR"]],
    cor.kIM = cor(discProps[["kIM"]], testProps[["kIM"]]),
    cor.MAR = cor(discProps[["MAR"]], testProps[["MAR"]])
  )
  if ("propVarExpl" %in% names(testProps)) {
    stats <- c(stats,
      propVarExpl = testProps[["propVarExpl"]],
      cor.kME = cor(discProps[["membership"]], testProps[["membership"]])
    )
  }
  stats
}

#' Network subset topology
#' 
#' Calculate the topological characteristics of a network subset.
#'  
#' @param adj Adjacency matrix for the network.
#' @param adjInd Indices of the network subset in \code{adj}.
#' @param dat (Optional) Underlying data for the network.
#' @param scaled (Optional) a row scaled \code{big.matrix} of \code{dat}.
#' @param datInd (Optional) Indices of the network subset in \code{dat}.
#'
#' @return
#'  A list of topological properties for the given network subset 
#' @seealso \code{\link[=subsetTestStats]{Between-network statistics}}
subsetProps <- function(adj, adjInd, dat=NULL, scaled=NULL, datInd=NULL) {
  props <- NetProps(adj@address, sort(adjInd))
  if (!is.null(dat)) {
    props <- c(props, DataProps(dat@address, scaled@address, datInd))
  }
  props
}
