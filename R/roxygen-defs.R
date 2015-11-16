#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param discCor a \code{\link[=bigMatrix-class]{bigMatrix}} containing the
#'   correlation structure among variables in the \emph{discovery} dataset.
#' @param testCor a \code{\link[=bigMatrix-class]{bigMatrix}} containing the
#'   correlation structure among variables in the \emph{test} dataset.
#' @param discIndices indices corresponding to the network module in 
#'  \code{discCor}.
#' @param testIndices indices corresponding to the network module in
#'  \code{testCor}.
#'
#' @name cor_params
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param adj a \code{\link[=bigMatrix-class]{bigMatrix}} containing an
#' adjacency matrix representation of the network (i.e. an \eqn{n * n} matrix 
#' containing the edge weights between each pair of variables).  
#' 
#' @name net_param
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param sge a \code{\link[=bigMatrix-class]{bigMatrix}} containing scaled 
#'  numeric data. Assumes columns correspond to variables (e.g. genes, 
#'  microbial operational taxonomic unit) and rows correspond to samples.
#'
#' @name dat_param
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
