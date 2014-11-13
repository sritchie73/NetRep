#' Fast wrapper functions for Rcpp functions
#' 
#' @template wrapper_desc
#' @template lowmem_inputs
#' @template no_sanity
#' @inheritParams coexp_params
#' @inheritParams adj_param
#' @inheritParams ge_param
#' @inheritParams ind_param
#' 
#' @name wrappers
NULL

#' @rdname wrappers
#' @return 
#'  \code{coexpStats:} a list containing the \emph{cor.coexp} and 
#'  \emph{mean.coexp} statistics for the specified module across datasets.
coexpStats <- function(
  discCoexp, discIndices, testCoexp, testIndices, lowmem=FALSE
) {
  if (lowmem) {
    discCoexp <- attach.big.matrix(discCoexp)
    testCoexp <- attach.big.matrix(testCoexp)
    on.exit({
      rm(discCoexp, testCoexp)
      gc()
    })
    poke(discCoexp, testCoexp)
  }
  res <- CoexpStats(
    discCoexp@address, discIndices, testCoexp@address, testIndices
  )
  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'  \code{netProps:} a list of network properties calculated from an adjacency
#'  matrix. These properties can either be scalers (summarising the whole
#'  network subset), or vectors (characterising some property for each node in
#'  the network subset).
adjProps <- function(adj, moduleIndices, lowmem=FALSE) {
  if (lowmem) {
    adj <- attach.big.matrix(adj)
    on.exit({
      rm(adj)
      gc()
    })
    poke(adj)
  }
  res <- AdjProps(adj@address, moduleIndices)
  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'   \code{dataProps:} a list of properties quantifying the relationship between
#'   a network subset and the underlying data the adjacency matrix was 
#'   calculated from. These properties can either be scalers (summarising the 
#'   whole network subset), or vectors (characterising some property for each 
#'   node in the network subset).
dataProps <- function(sge, moduleIndices, lowmem=FALSE) {
  if (lowmem) {
    sge <- attach.big.matrix(sge)
    on.exit({
      rm(sge)
      gc()
    })
    poke(sge)
  } 
  res <- DataProps(sge@address, moduleIndices)
  lapply(res, as.vector)
}
  
