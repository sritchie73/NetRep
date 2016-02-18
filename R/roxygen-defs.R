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
#'  to in the discovery dataset. If there are multiple discovery datasets 
#'  then this argument should be a list of such vectors.
#' @param modules a list of vectors, one for each \code{discovery} dataset, 
#'  of modules to perform the analysis on. The default is to analyse all modules
#'  with the exception of those specified in \code{backgroundLabel}.
#' @param backgroundLabel a single label that nodes that do not belong to any
#'  module are assigned. The default is "0".
#' @param discovery a vector of names or indices denoting the discovery dataset(s).
#' @param test a list of vectors of names or indices denoting the test datasets
#'  for each \code{discovery} dataset. Alternatively can be provided as vector
#'  if the test datasets are the same for all 'discovery' datasets (e.g. for 
#'  performing a pairwise comparison).
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#' 
#' @name common_params
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param simplify logical; if \code{TRUE}, simplify the structure of the output
#'  list if possible (see Return Value).
#'  
#' @name simplify_param
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param nCores number of cores to parallelise the calculation of network 
#'  properties over. Ignored if the user has already registered a parallel 
#'  backend.
#'
#' @name par_param
NULL