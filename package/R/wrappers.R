#' Fast wrapper functions for Rcpp functions
#' 
#' @template wrapper_desc
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
#'  \code{netProps:} a list of network properties calculated from an adjacency
#'  matrix. These properties can either be scalers (summarising the whole
#'  network subset), or vectors (characterising some property for each node in
#'  the network subset).
adjProps <- function(adj, moduleIndices) {
  # Attach the big.matrix object if not attached yet
  is.attached <- adj@attached
  if (!is.attached)
    adj <- attach.bigMatrix(adj)
  
  res <- AdjProps(adj@matrix@address, moduleIndices)
  
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
