#' @param geneExpression the gene expression matrix for the dataset of interest,
#'   or a list of gene expression matrices, one for each dataset. Columns are 
#'   expected to be genes, rows samples. See details for expected input types.
#' @param coexpression the coexpression matrix for the dataset of interest, or a
#'   list of coexpression matrices, one for each dataset. See details for
#'   expected input types.
#' @param adjacency the adjacecny matrix for the dataset of interest, or a
#'   list of adjacency matrices, one for each dataset. See details for
#'   expected input types.
#' @param moduleAssignments a vector assigning genes to modules, or a list of 
#'  such vectors. See details.
#' @param modules a vector of modules to apply the function to. See details.
#' @param discovery name or index denoting which dataset the module of
#'  interest was discovered in. See details.
#' @param test name or index denoting which dataset to apply the function to.
#'  See details.
#' @param ... additional arguments used when reading in a matrix from file. See
#'  \code{\link{read.bigMatrix}}.
#' 
#' @details
#'  The matrices containing the \code{geneExpression}, \code{coexpression}, and
#'  \code{adjacency} should be \code{\link{bigMatrix}} object(s), but the
#'  function will also work with regular matrices, files containing the matrix
#'  data, or file paths to the descriptor files for a 'bigMatrix', in which 
#'  case temporary 'bigMatrix' objects will be created. It is advisable to keep
#'  large matrices stored as 'bigMatrix' objects, as the conversion process 
#'  involves writing out the matrix as binary data to disk. 'bigMatrix' objects
#'  have the added advantage of having instant load times in new R sessions, and
#'  can be accessed and used the same way as regular matrices. See
#'  \code{\link{bigMatrix}} for more details.
#'  
#'  This function can be used in three ways. First, the \code{geneExpression}, 
#'  \code{coexpression}, and \code{adjacency} can be provided as lists of 
#'  matrices, where each element corresponds to a dataset of interest. In this
#'  case, the function will be applied to the \code{test} dataset for the subset
#'  of genes corresponding to modules discovered in the \code{discovery} dataset.
#'  By default, the function is applied on dataset 1, on the modules discovered 
#'  in dataset 1. Secondly, the function also accepts single matrices for each,
#'  simplifying the input for users with only one dataset of interest. Finally,
#'  \code{modules} and \code{moduleAssignments} may be omitted, in which case
#'  the function is applied on all genes present in the \code{geneExpression}, 
#'  \code{coexpression}, and \code{adjacency}. This may be useful when 
#'  subsetting each matrix by a custom gene set of interest.
#' 
#' @name api_inputs
NULL
