#' Fast permutation procedure for testing network module replication
#'
#' Functions for assessing the replication/preservation of network
#' topology for weighted gene coexpression network modules in one or more
#' independent datasets through permutation testing.
#' 
#' The main function for this package is \code{\link{modulePreservation}}. Other
#' useful functions include \code{\link{networkProperties}} for calculating the
#' topological properties of a module, and \code{\link{plotModule}} and other
#' plotting functions for visualising a module. 
#'
#' @docType package
#' @name netrep
#' @useDynLib netrep
#' @importFrom Rcpp evalCpp
#' @import BH
#' @import bigmemory
NULL
