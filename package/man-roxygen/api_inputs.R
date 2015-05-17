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
#' If any of the \code{geneExpression}, \code{coexpression}, \code{adjacency}, 
#' or \code{moduleAssignments} arguments are list objects, then they all must 
#' be. In this case, each list element corresponds to a dataset that the 
#' respective network components have been calculated in. Providing the 
#' \code{geneExpression} is optional, but limits the scope of the analysis or 
#' function output. The \code{discovery} argument denotes which dataset the 
#' \code{modules} of interest were discovered in, and the \code{test} argument 
#' controls which dataset the function is calculated on. Alternatively, the
#' \code{geneExpression}, \code{coexpression}, \code{adjacency}, and
#' \code{moduleAssignments} can be provided only for the discovery dataset.
#' 
#' @name api_inputs
NULL
