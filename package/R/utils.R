#' Value Matching and Subsetting
#' 
#' @description
#' This set of functions provides shortcuts for value matching and subsetting,
#' on top of the functionality provided by \code{\link[base]{\%in\%}}.
#' 
#' \code{\%nin\%} returns a logical vector indicating if elements of \code{x} 
#' are not in \code{table}, This is the opposite of \code{\%in\%}.
#' 
#' \code{\%sub_in\%} returns the elements \code{x} that are \code{\%in\%} 
#' \code{table} rather than a logical vector.
#' 
#' \code{\%sub_nin\%} returns the elements \code{x} that are \code{\%nin\%} 
#' \code{table} rather than a logical vector.
#' 
#' @param x vector or \code{NULL}: the values to be matched. \link[base]{Long vectors} 
#'   are supported.
#' @param table vector or \code{NULL}: the values to be matched against. 
#'   \link[base]{Long vectors} are not supported.
#' @name matchsub
NULL

#' @rdname matchsub
#' @export
`%nin%` <- function(x, table) !(x %in% table)

#' @rdname matchsub
#' @export
`%sub_in%` <- function(x, table) x[x %in% table]

#' @rdname matchsub
#' @export
`%sub_nin%` <- function(x, table) x[x %nin% table]

#' Verbose Concatenate and Print
#' 
#' Concatenate and output the objects only if the \code{verbose} flag is set 
#' to \code{TRUE}.
#' 
#' @param verbose logical. If \code{TRUE}, passes the rest of the arguments to
#'   \code{\link{cat}}
#' @param ... arguments to pass to \code{\link{cat}}
vCat <- function(verbose, ...) if(verbose) cat(..., "\n")

#' Binds two dimensional arrays along the third dimension.
#' @importFrom abind abind
abind3 <- function(...) abind(..., along=3)
