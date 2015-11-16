#' Fast permutation procedure for testing network module replication
#'
#' Functions for assessing the replication/preservation of a network module's
#' topology across datasets through permutation testing. This method is suitable
#' for interaction networks that are not sparse. These include gene coexpression
#' networks, protein-protein interaction networks, and microbial interaction 
#' networks. Modules may be defined as tightly connected clusters or through 
#' ontologies. Application of this method can answer questions such as; (1) does
#' a disease-associated gene coexpression module replicate in an independent 
#' cohort? (2) are these gene coexpression modules preserved across tissues or 
#' tissue specific? (3) are these modules conserved across species? (4) are 
#' microbial communities preseved across multiple locations?
#' 
#' The main function for this package is \code{\link{modulePreservation}}. 
#' Several functions for downstream are also provided: 
#' \code{\link{networkProperties}} for calculating the topological properties 
#' of a module, and \code{\link{plotModule}} for visualising a module.
#'
#' @docType package
#' @name NetRep
#' @useDynLib NetRep
#' @importFrom Rcpp evalCpp
#' @import BH
#' @import bigmemory
NULL
