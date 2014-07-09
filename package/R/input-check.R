#' Helper function for checking input sanity
#' 
#' Throws an error if \code{netSets} is not ok.
#' 
#' @param netSets a list of networks
checkNetSets <- function(netSets) {
  stopifnot(sapply(netSets, class) %in% c("big.matrix", "NULL"))
  stopifnot(sapply(netSets[!is.null(netSets)], function(x) !is.null(rownames(x))))
  stopifnot(sapply(netSets[!is.null(netSets)], function(x) !is.null(colnames(x))))
}

#' Helper function for checking input sanity
#' 
#' Throws an error if \code{datSets} is not ok.
#' 
#' @param datSets a list of data matrices
checkDatSets <- function(datSets) {
  stopifnot(sapply(datSets, class) %in% c("big.matrix", "NULL"))
  stopifnot(sapply(datSets[!is.null(datSets)], function(x) !is.null(rownames(x))))
}

#' Check consistency of all dataset lists 
#' 
#' Throws an error if the dataset lists are not consistent.
#' 
#' @param datSets a list of data matrices (or \code{NULL})
#' @param netSets a list of networks (or \code{NULL})
checkSetConsistency <- function(datSets, netSets) {
  
}