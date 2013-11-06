#' Run permutations on the provided function
#' 
#' Does \code{nPerm} calculations on the provided expression in parallel. The 
#' \code{permuteVector} is shuffled, and passed to the provided function.
#' 
#' @param FUN the function to pass the \code{permuteVector} to.
#' @param permuteVector vector to shuffle. Names of the vector will remain
#'  constant for each permutation.
#' @param nPerm Number of permutations to run
#' @param verbose Logical; If TRUE, writes the current progress to file at 
#'  every iteration.
#' @param logFile The file to output progress to if \code{verbose} is 
#'  \code{TRUE}.
#' @param ... optional arguments to \code{\link{foreach}}
#' @export
doParPermute <- function(FUN, permuteVector, nPerm=1000, 
                         verbose=TRUE, logFile='log.txt', ...) {
  require(foreach)
  original <- permuteVector
  foreach (i = 1:nPerm, ...) %dopar% {
    permuteVector <- sample(original, length(original))
    if (!is.null(names(original))) {
      names(permuteVector) <- names(original)
    }
    if (verbose) {
      LogProgress(i, nPerm, file=logFile)
    }
    FUN(permuteVector)
  }
}

