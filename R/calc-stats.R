#' Calculate module preservation statistics.
#' 
#' Calculates the modules preservation statistics using the pre-calculated 
#' module properties along with the correlation structure matrices.
#' 
#' @param discProps module properties in the \emph{discovery} dataset calculated 
#'  using \code{\link{netProps}} and \code{\link{dataProps}}.
#' @param testProps module properties in the \emph{test} dataset calculated 
#'  using \code{\link{netProps}} and \code{\link{dataProps}}.
#' @param discCor 'bigMatrix' containing the correlation structure matrix 
#'  for the \emph{discovery} dataset.
#' @param discIndices module indices in the discovery dataset.
#' @param testCor 'bigMatrix' containing the correlation structure matrix 
#'  for the \emph{test} dataset.
#' @param testIndices module indices in the test dataset.
#'  
#' @section Design decision:
#'  Most of the module presevation statistics are calculated as the correlation 
#'  of some network property between the \emph{discovery} and \emph{test} 
#'  datasets. The permutation procedure randomly samples from the \emph{test} 
#'  dataset, so it therefore makes sense to calculate the network properties in 
#'  the \emph{discovery} dataset once only. The code is thus split into 
#'  calculations of the network properties, and then subsequent calculation of 
#'  the statistics. The exception to this design are the correlation structure 
#'  based statistics as they require all entries of the correlation structure 
#'  matrices for that module. As modules can be quite large (e.g., 10,000 
#'  variables), this would be dramatically increase memory overhead to apply 
#'  this strategy.
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
#' statistics: \emph{density}, \emph{cor.kIM}, \emph{pve}, \emph{mean.MM},
#' \emph{cor.MM}, \emph{cor.cor}, and \emph{mean.cor}.
calcStats <- function(
  discProps, testProps, discCor, discIndices, testCor, testIndices
) {
  cList <- corStats(discCor, discIndices, testCor, testIndices)
  stats <- c(
    density = testProps[["density"]],
    cor.kIM = cor(discProps[["kIM"]], testProps[["kIM"]]),
    cor.cor = cor(cList[["cor.discovery"]], cList[["cor.test"]]),
    mean.cor = cList[["mean.cor"]]
  )
  if ("pve" %in% names(testProps)) { # Detect if data has been provided
    stats <- c(
      stats, propVarExpl = testProps[["pve"]],
      mean.MM = mean(sign(discProps[["MM"]]) * testProps[["MM"]]),
      cor.MM = cor(discProps[["MM"]], testProps[["MM"]])
    )
  }
  stats
}
