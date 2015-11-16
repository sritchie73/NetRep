#' Fast wrapper functions for Rcpp functions
#' 
#' Wrapper functions for interfacing with the Rcpp implementations of 
#' \code{\link{CoexpStats}}, \code{\link{AdjProps}}, and \code{\link{DataProps}}.
#' They provide a nicer interface for input handling and format the output in
#' an R-friendly way.
#'  
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
  discCoexp, discIndices, testCoexp, testIndices
) {
  # Attach the big.matrix objects if not attached yet
  disc.attached <- discCoexp@attached
  test.attached <- testCoexp@attached
  if (!disc.attached)
    discCoexp <- attach.bigMatrix(discCoexp)
  if (!test.attached)
    testCoexp <- attach.bigMatrix(testCoexp)
  
  res <- CoexpStats(
    discCoexp@matrix@address, discIndices, 
    testCoexp@matrix@address, testIndices
  )
  
  # detach big.matrix objects if they were detached to begin with
  if (!disc.attached)
    discCoexp <- detach.bigMatrix(discCoexp)
  if (!test.attached)
    testCoexp <- detach.bigMatrix(testCoexp)

  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'  \code{netProps:} a list of containing:  
#'  \enumerate{
#'    \item{kIM:}{
#'      The intramodular connectivity (the weighted within-subset 
#'      degree for each network node)
#'    }
#'    \item{density:}{
#'      The mean edge weight within the network module 
#'    }
#'  }
netProps <- function(adj, moduleIndices) {
  # Attach the big.matrix object if not attached yet
  is.attached <- adj@attached
  if (!is.attached)
    adj <- attach.bigMatrix(adj)
  
  res <- NetProps(adj@matrix@address, moduleIndices)
  
  # detach big.matrix objects if they were detached to begin with
  if (!is.attached)
    adj <- detach.bigMatrix(adj)
  
  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'   \code{dataProps:} a list of properties quantifying the relationship between
#'   a network subset and the underlying data the adjacency matrix was 
#'   calculated from. These properties can either be scalers (summarising the 
#'   whole network subset), or vectors (characterising some property for each 
#'   node in the network subset).
dataProps <- function(sge, moduleIndices) {
  # Attach the big.matrix object if not attached yet
  is.attached <- sge@attached
  if (!is.attached)
    sge <- attach.bigMatrix(sge)
  
  res <- DataProps(sge@matrix@address, moduleIndices)
  
  # detach big.matrix objects if they were detached to begin with
  if (!is.attached)
    sge <- detach.bigMatrix(sge)
  
  lapply(res, as.vector)
}

#' @rdname wrappers
#' 
#' @description
#'  combines \code{dataProps} and \code{adjProps}
moduleProps <- function(adj, moduleIndices, sge) {
  geProps <- NULL
  if (!is.null(sge))
    geProps <- dataProps(sge, moduleIndices)
  adjProps <- adjProps(adj, moduleIndices)
  
  c(geProps, adjProps)
}
