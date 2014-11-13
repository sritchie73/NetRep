#' Calculate module preservation statistics
#' 
#' @description
#'  A set of functions for calculating the module preservation statistics, and 
#'  their underlying network properties, for a single module.
#'  
#' @section Design Rationale:
#'  Most of the module presevation statistics are calculated as the correlation 
#'  of some network properties across the \emph{discovery} and \emph{test}
#'  datasets. The permutation procedure randomly samples from the \emph{test}
#'  dataset, so it therefore makes sense to calculate the network properties in
#'  the \emph{discovery} dataset once. 
#'  
#'  The code is thus split into calculations of the network properties, and then
#'  calculating the statistics from these. The ability to calculate the network
#'  properties is also useful for downstream analysis.
#'  
#'  The exceptions to this design are the coexpression based statistics. The 
#'  \emph{cor.coexp} and \emph{mean.coexp} require all entries of the pairwise
#'  coexpression matrix for the module. As modules can be quite large (e.g.,
#'  10,000 probes), this would be prohibitive to pass around in memory.
#' 
#' @template Langfelder-ref
#' 
#' @template no_sanity
#' @template lowmem_inputs
#' 
#' @name module-preservation-stats
#' 
NULL

#' @rdname module-preservation-stats
#' 
#' @description
#'  \code{moduleProps}: calculate a module's network properties. 
#'  
#' @details
#'  \code{moduleProps}: providing gene expression data is optional. If not 
#'  provided, only the \emph{kIM} and \emph{mean.adj} are calculated.
#'  
#' @inheritParams ge_param
#' @inheritParams ind_param
#' @inheritParams adj_param
#' 
#' @return
#'  \code{moduleProps}: A list containing the following network properties:
#'  \enumerate{
#'    \item{\emph{"kIM"}}:{
#'      A vector whose entries contain the intramodular connectivity for each
#'      gene in the module.
#'    }
#'    \item{\emph{"mean.adj"}}:{
#'      The density of the module: the average pairwise-gene adjacency.
#'    }
#'    \item{\emph{"SEP"}}:{
#'      A numeric vector denoting the summary expression profile for the module.
#'    }
#'    \item{\emph{"MM"}}:{
#'      A vector whose entries contain the module membership for each gene; the
#'      correlation between each gene and the summary expression profile.
#'    }
#'    \item{\emph{"pve"}}:{
#'      The proportion of variance explained in the module's underlying gene 
#'      expression by the summary expression profile.
#'    }
#'  }
#'  
moduleProps <- function(adj, moduleIndices, sge=NULL, lowmem=FALSE) {
  props <- adjProps(adj, moduleIndices, lowmem)
  if (!is.null(sge)) {
    props <- c(props, dataProps(sge, moduleIndices, lowmem))
  }
  props
}

#' @rdname module-preservation-stats
#' 
#' @description
#'  \code{calcSplitTestStats}: calculate the module preservation statistics 
#'  from two lists of network properties obtained from \code{moduleProps}.
#' 
#' @param discProps properties of the network module in the \emph{discovery}
#'  dataset (calculated by \code{moduleProps}).
#' @param testProps properties of the network module in the \emph{test} dataset
#'  (calculated by \code{moduleProps}). 
#' 
#' @return 
#'  \code{calcSplitTestStats}: A vector containing the \emph{mean.adj}, and
#'  \emph{cor.kIM}, as well as the \code{pve}, \code{mean.MM}, and \code{cor.MM}
#'  if gene expression data has been supplied when calculating \code{testProps} 
#'  and \code{discProps}.
calcSplitTestStats <- function(discProps, testProps) {
  stopifnot(is.list(discProps) & is.list(testProps))
  stopifnot(length(discProps) == length(testProps))
  
  stats <- c(
    mean.adj = testProps[["mean.adj"]],
    cor.kIM = cor(discProps[["kIM"]], testProps[["kIM"]])
  )
  if ("pve" %in% names(testProps)) { # Detect if data has been provided
    stats <- c(stats,
      pve = testProps[["pve"]],
      mean.MM = mean(sign(discProps[["MM"]]) * testProps[["MM"]]),
      cor.MM = cor(discProps[["MM"]], testProps[["MM"]])
    )
  }
  stats
}


#' @rdname module-preservation-stats
#' 
#' @description
#'   \code{calcSharedTestStats}: calculate remaining module preservation 
#'   statistics for which network properties have not been calculated in 
#'   advance.
#'
#' @inheritParams coexp_params
#' 
#' @return 
#'  \code{calcSharedTestStats}: A vector of containing the \emph{cor.coexp} and
#'  \emph{mean.coexp} statistics for the specified module.
calcSharedTestStats <- function(
  discCoexp, discIndices, testCoexp, testIndices, lowmem
) {
  unlist(coexpStats(discCoexp, discIndices, testCoexp, testIndices, lowmem))
}

#' @rdname module-preservation-stats
#' @description
#'  \code{calcStats}: interface function that calls \code{calcSplitTestStats}
#'  and \code{calcSharedTestStats} and combines their results into a single
#'  vector.
#'  
#' @return
#'  \code{calcStats}: a vector containing the values for all the module 
#'  preservation statistics: \emph{mean.adj}, \emph{cor.kIM}, \emph{pve},
#'  \emph{mean.MM}, \emph{cor.MM}, \emph{cor.coexp}, and \emph{mean.coexp}.
calcStats <- function(
  discProps, testProps, discCoexp, discIndices, testCoexp, testIndices, lowmem
) {
  c(
    calcSplitTestStats(discProps, testProps),
    calcSharedTestStats(
      discCoexp, discIndices, testCoexp, testIndices, lowmem
    )
  )
}
