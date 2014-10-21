#' Scale and Center the rows of a 'big.matrix'
#' 
#' Create a new scaled \code{\link[bigmemory]{big.matrix}}.  
#' 
#' @param x the \code{\link[bigmemory]{big.matrix}} to scale.
#' @param backingfile,backingpath,descriptorfile backingfile information for the
#'   new scaled big.matrix to be stored in, see 
#'   \code{\link[bigmemory]{big.matrix}}.
#' @return
#'  A new \code{\link[bigmemory]{big.matrix}} containing the scaled and centered
#'  rows of \code{x}.
#' 
#' @seealso \code{\link[base]{scale}} \code{\link[bigmemory]{big.matrix}}
#' @note
#'  Unlike \code{\link[base]{scale}} this function applies the scale function to
#'  the rows of the matrix instead of the columns.
#'  
#'  This function has only been implemented for 
#'  \code{\link[bigmemory]{big.matrix}} objects of \code{type = "double"}. 
#'  
#'  To apply the scale function in place without creating a new results matrix,
#'  simply call \code{Scale(x@@address, x@@address)}.
#' 
scaleBigMatrix <- function(
  x, backingfile=NULL, backingpath=NULL, descriptorfile=NULL
) {
  res <- big.matrix(nrow(x), ncol(x), typeof(x), NULL, dimnames(x), FALSE,
                    backingfile, backingpath, descriptorfile)
  Scale(x@address, res@address)
  res
}

#' Get the range of a big.matrix or its subset
#' 
#' Note subsetting only applies to columns: this is only meant for use with the 
#' gene expression matrix.
#' 
#' @param x a big matrix
#' @param subsetIndices an optional vector to subset the matrix by
#' 
rangeBigMatrix <- function(x, subsetIndices) {
  if (missing(subsetIndices)) {
    res <- BigRange(x@address)
  } else {
    res <- RangeSubset(x@address, subsetIndices)
  }
  unlist(lapply(res, as.vector))
}

#' Check if all entries of a `big.matrix` are Finite
#' 
#' If there are non-finite entires (\code{NA}, \code{NaN}, \code{-Inf}, 
#' \code{Inf}), throw an exception.
#' 
#' @param x a \code{\link[bigmemory]{big.matrix}}
checkFinite <- function(x) {
  CheckFinite(x@address)
}
