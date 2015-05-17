#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param discCoexp a matrix (see Inputs section) of the pairwise coexpression 
#'   for the \emph{discovery} dataset.
#' @param testCoexp a matrix (see Inputs section) of the pairwise coexpression 
#'   for the \emph{test} dataset.
#' @param discIndices indices corresponding to the network module in the matrix
#'   of pairwise-coexpression in the \emph{discovery} dataset.
#' @param testIndices indices corresponding to the network module, or a random 
#'   sample of the same size, in the matrix of pairwise-coexpression in the 
#'   \emph{test} dataset.
#'   
#' @name coexp_params
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param adj a matrix (see Inputs section) of pairwise gene adjacencies.
#' 
#' @name adj_param
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param sge a matrix (see Inputs section) of gene expression data. Assumes 
#'  columns are genes/probes, rows are samples, and that probes have been scaled 
#'  (This can be done using \code{\link{scaleBigMatrix}}).
#'
#' @name ge_param
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param moduleIndices indices for the network module of interest in the
#'   supplied data matrix.
#'
#' @name ind_param
NULL
