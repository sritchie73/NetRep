#' Scale and Center the rows of a 'big.matrix'
#' 
#' Create a new scaled \code{\link[bigmemory]{big.matrix}}.  
#' 
#' @param x the matrix to scale.
#' @param obj.dir the directory to store the object in
#' @return
#'  A descriptor to the scaled version of x.
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
scaleBigMatrix <- function(x, obj.dir) {
  bigx <- dynamicMatLoad(x)
  on.exit({
    rm(bigx, res)
    gc()  
  })
  poke(bigx)
  

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
    nrow(bigx), ncol(bigx), typeof(bigx), NULL, dimnames(x), FALSE,
    backingfile, obj.dir, descriptorfile
  )
  Scale(bigx@address, res@address)

  return(file.path(obj.dir, descriptorfile))
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
  bigx <- dynamicMatLoad(x)
  on.exit({
    rm(bigx)
    gc()
  })
  poke(bigx)
  CheckFinite(bigx@address)
  return(NULL)
}
