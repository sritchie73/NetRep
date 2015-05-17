#' @param geneExpression a \code{\link{bigMatrix}} object, or a list of 
#'  \code{bigMatrix} objects, see details.
#' @param coexpression a \code{\link{bigMatrix}} object, or a list of 
#'  \code{bigMatrix} objects, see details.
#' @param adjacency a \code{\link{bigMatrix}} object, or a list of 
#'  \code{bigMatrix} objects, see details.
#' @param moduleAssignments a vector assigning genes to modules, or a list of 
#'  such vectors. See details.
#' @param modules a vector of modules to apply the function to. See details.
#' @param discovery name or index denoting which dataset the module of
#'  interest was discovered in. See details.
#' @param test name or index denoting which dataset to apply the function to.
#'  See details.
#' 
#' @details
#' If any of the arguments \code{geneExpression}, \code{coexpression}, 
#' \code{adjacency}, or \code{modules} are provided as 
#' 
#' @return
#' A list of network properties 
#' 
#' @name api_inputs
NULL