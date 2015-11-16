#' Fast wrapper functions for Rcpp functions
#' 
#' Wrapper functions for interfacing with the Rcpp implementations of 
#' \code{\link{CorStats}}, \code{\link{NetProps}}, and \code{\link{DataProps}}.
#' They provide a nicer interface for input handling and format the output in
#' an R-friendly way.
#'  
#' @inheritParams cor_params
#' @inheritParams net_param
#' @inheritParams dat_param
#' @inheritParams ind_param
#'
#' @references
#'  \enumerate{
#'     \item{
#'       Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my
#'       network module preserved and reproducible?} PLoS Comput. Biol. 
#'       \strong{7}, e1001057 (2011). 
#'     }
#'  }
#'
#' @name wrappers
NULL

#' @rdname wrappers
#' @return 
#'  \code{corStats:} a list containing:
#'  \enumerate{
#'    \item{\emph{cor.discovery}:}{
#'       A flattened vector of the module's correlation structure in the 
#'       \emph{discovery} dataset.
#'    } \item{\emph{cor.test}:}{
#'       A flattened vector of the module's correlation structure in the 
#'       \emph{test} dataset.
#'    } \item{\emph{mean.cor}:}{
#'      The mean sign-aware correlation density of the network subset.
#'    }
#'  } 
corStats <- function(
  discCor, discIndices, testCor, testIndices
) {
  # Attach the bigMatrix objects if not attached yet
  disc.attached <- discCor@attached
  test.attached <- testCor@attached
  if (!disc.attached)
    discCor <- attach.bigMatrix(discCor)
  if (!test.attached)
    testCor <- attach.bigMatrix(testCor)
  
  res <- CorStats(
    discCor@matrix@address, discIndices, 
    testCor@matrix@address, testIndices
  )
  
  # detach bigMatrix objects if they were detached to begin with
  if (!disc.attached)
    discCor <- detach.bigMatrix(discCor)
  if (!test.attached)
    testCor <- detach.bigMatrix(testCor)

  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'  \code{netProps:} a list of containing:  
#'  \enumerate{
#'    \item{\emph{kIM}:}{
#'      The intramodular connectivity (the weighted within-subset 
#'      degree for each network node)
#'    }
#'    \item{\emph{density}:}{
#'      The mean edge weight within the network module 
#'    }
#'  }
netProps <- function(adj, moduleIndices) {
  # Attach the bigMatrix object if not attached yet
  is.attached <- adj@attached
  if (!is.attached)
    adj <- attach.bigMatrix(adj)
  
  res <- NetProps(adj@matrix@address, moduleIndices)
  
  # detach bigMatrix objects if they were detached to begin with
  if (!is.attached)
    adj <- detach.bigMatrix(adj)
  
  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'   \code{dataProps:} a list containing:
#'   \enumerate{
#'    \item{\emph{moduleSummary}:}{The module summary vector.}
#'    \item{\emph{MM}:}{
#'      The module membership for each variable, i.e. the correlation between
#'      each variable and the module summary vector.
#'    }
#'    \item{\emph{pve}:}{
#'      The proportion of module variance explained by the module summary vector. 
#'    }
#'   }
dataProps <- function(sdat, moduleIndices) {
  # Attach the bigMatrix object if not attached yet
  is.attached <- sdat@attached
  if (!is.attached)
    sdat <- attach.bigMatrix(sdat)
  
  res <- DataProps(sdat@matrix@address, moduleIndices)
  
  # detach bigMatrix objects if they were detached to begin with
  if (!is.attached)
    sdat <- detach.bigMatrix(sdat)
  
  lapply(res, as.vector)
}

#' @rdname wrappers
#' 
#' @description
#'  combines \code{dataProps} and \code{netProps}
moduleProps <- function(adj, moduleIndices, sdat) {
  datProps <- NULL
  if (!is.null(sdat))
    datProps <- dataProps(sdat, moduleIndices)
  netProps <- netProps(adj, moduleIndices)
  
  c(datProps, netProps)
}
