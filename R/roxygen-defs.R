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

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param data a list of numeric matrices. Each entry of the list corresponds to
#'   a dataset and contains the data used to infer the interaction network
#'   between variables (e.g. genes). Expects matrix columns to correspond to
#'   variables and matrix rows to correspond to samples.
#' @param correlation a list of matrices. Each entry of the list corresponds to a 
#'  dataset and contains an \eqn{n * n} matrix of the correlation between 
#'  each pair of variables in the dataset.
#' @param network a list of matrices. Each entry of the list corresponds to a 
#'  dataset and contains an \eqn{n * n} matrix of the network edge weights 
#'  between each pair of variables in the dataset.
#' @param moduleAssignments a vector containing the module each variable belongs
#' to in the discovery dataset. If there are multiple discovery datasets 
#' then this argument should be a list of such vectors.
#' @param discovery name or index denoting which the discovery dataset.
#' @param test name or index denoting the dataset to analyse the module in.
#'
#' @name common_params
NULL

