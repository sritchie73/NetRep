#' Calculate module preservation statistics
#' 
#' Calculates the modules preservation statistics using the sets of 
#' pre-calculated properties, along with the coexpression networks.
#' 
#' @param discProps network properties calculated for the module from the
#'  discovery dataset using \code{\link{adjProps}} and \code{\link{dataProps}}.
#' @param testProps network properties calculated for the module from the
#'  test dataset using \code{\link{adjProps}} and \code{\link{dataProps}}.
#' @param discCoexp 'bigMatrix' containing the coexpression network for the discovery dataset.
#' @param discIndices module indices in the discovery dataset.
#' @param testCoexp 'bigMatrix' containing the coexpression network for the test dataset.
#' @param testIndices module indices in the test dataset.
#'  
#' @section Design decision:
#'  Most of the module presevation statistics are calculated as the correlation 
#'  of some network properties across the \emph{discovery} and \emph{test} 
#'  datasets. The permutation procedure randomly samples from the \emph{test} 
#'  dataset, so it therefore makes sense to calculate the network properties in 
#'  the \emph{discovery} dataset once only. The code is thus split into 
#'  calculations of the network properties, and then calculating the statistics 
#'  from these. The exceptions to this design are the coexpression based
#'  statistics. The \emph{cor.coexp} and \emph{mean.coexp} require all entries
#'  of the pairwise coexpression matrix for the module. As modules can be quite
#'  large (e.g., 10,000 probes), this would be prohibitive to pass around in
#'  memory.
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
#' @return 
#' a vector containing the values for all the module preservation
#' statistics: \emph{mean.adj}, \emph{cor.kIM}, \emph{pve}, \emph{mean.MM},
#' \emph{cor.MM}, \emph{cor.coexp}, and \emph{mean.coexp}.
calcStats <- function(
  discProps, testProps, discCoexp, discIndices, testCoexp, testIndices
) {
  stats <- c(
    mean.adj = testProps[["mean.adj"]],
    cor.kIM = cor(discProps[["kIM"]], testProps[["kIM"]]),
    unlist(coexpStats(discCoexp, discIndices, testCoexp, testIndices))
  )
  if ("pve" %in% names(testProps)) { # Detect if data has been provided
    stats <- c(
      stats, pve = testProps[["pve"]],
      mean.MM = mean(sign(discProps[["MM"]]) * testProps[["MM"]]),
      cor.MM = cor(discProps[["MM"]], testProps[["MM"]])
    )
  }
  stats
}
