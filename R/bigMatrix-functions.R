#' Scale and Center the rows of a 'big.matrix'
#' 
#' Create a new \code{\link[=bigMatrix-class]{bigMatrix}} containing scaled gene
#' expression. Genes are expected to be columns.
#' 
#' @param x the bigMatrix to scale.
#' @param tmp.dir temporary directory to store the backingfiles for the 
#'  bigMatrix.
#' @return
#'  an object of type \code{bigMatrix}.
#' 
#' @seealso \code{\link[base]{scale}} \code{\link[bigmemory]{big.matrix}}
#' @note
#'  The caller is responsible for cleaning up the temporary directory.
#' 
scaleBigMatrix <- function(x, tmp.dir) {
  is.attached <- x@attached
  if (!is.attached)
    x <- attach.bigMatrix(x)  

  while (TRUE) {
    stamp <- getUUID()
    descriptorfile <- paste0("scaled", stamp, ".desc")
    backingfile <- paste0("scaled", stamp, ".bin")
    
    # Handle the infintesimally small chance of a UUID collision
    if (!file.exists(descriptorfile) & !file.exists(backingfile)) {
      break 
    }
  }
  
  res <- big.matrix(
    nrow(x), ncol(x), typeof(x), NULL, dimnames(x), FALSE,
    backingfile, tmp.dir, descriptorfile
  )
  Scale(x@matrix@address, res@address)
  
  if (is.attached)
    x <- detach.bigMatrix(x)

  suppressWarnings(
    load.bigMatrix(file.path(tmp.dir, gsub(".desc", "", descriptorfile)))
  )  
}

#' Get the range of a bigMatrix or its subset
#' 
#' Note subsetting only applies to columns: this is only meant for use with the 
#' gene expression matrix.
#' 
#' @param x a bigMatrix
#' @param subsetIndices an optional vector to subset the matrix by
#' 
rangeBigMatrix <- function(x, subsetIndices) {
  is.attached <- x@attached
  if (!is.attached)
    x <- attach.bigMatrix(x)
  
  if (missing(subsetIndices)) {
    res <- BigRange(x@matrix@address)
  } else {
    res <- RangeSubset(x@matrix@address, subsetIndices)
  }
  
  if (is.attached)
    x <- detach.bigMatrix(x)
  
  unlist(lapply(res, as.vector))
}

#' Check if all entries of a `bigMatrix` are Finite
#' 
#' If there are non-finite entires (\code{NA}, \code{NaN}, \code{-Inf}, 
#' \code{Inf}), throw an exception. 
#' 
#' @param x a \code{\link[=bigMatrix-class]{bigMatrix}}
checkFinite <- function(x) {
  is.attached <- x@attached
  if (!is.attached)
    x <- attach.bigMatrix(x)
  
  tryCatch({
    CheckFinite(x@matrix@address)
  }, error=function(e) {
    stop(e$message)
  }, finally = {
    if (is.attached)
      x <- detach.bigMatrix(x)
  })
  
  return(NULL)
}
