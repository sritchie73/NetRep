#' Check if all entries of a 'bigMatrix' are finite
#' 
#' If there are non-finite entires (\code{NA}, \code{NaN}, \code{-Inf}, 
#' \code{Inf}), throw an exception. 
#' 
#' @param x a \code{\link{bigMatrix}}
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
