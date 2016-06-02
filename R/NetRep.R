#' Fast permutation procedure for testing network module replication
#' 
#' Functions for assessing the replication/preservation of a network module's 
#' topology across datasets through permutation testing. This is suitable for 
#' networks that can be meaningfully inferred from multiple datasets. These 
#' include gene coexpression networks, protein-protein interaction networks, and
#' microbial interaction networks. Modules within these networks consist of 
#' groups of nodes that are particularly interesting: for example a group of 
#' tightly connected genes associated with a disease, groups of genes annotated
#' with the same term in the Gene Ontology database, or groups of interacting
#' microbial species, i.e. communities. Application of this method can answer
#' questions such as; (1) do the relationships between genes in a module 
#' replicate in an independent cohort? (2) are these gene coexpression modules
#' preserved across tissues or tissue specific? (3) are these modules conserved
#' across species? (4) are microbial communities preseved across multiple spatial
#' locations?
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
#' @import RcppArmadillo
#' @import BH
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("NetRep", libpath)
}
