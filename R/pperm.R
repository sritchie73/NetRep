#' @title Permutation Test
#' 
#' @description 
#'  Evaluates the statistical significance of a test statistic given a vector
#'  of "nulls": values for that test statistic drawn from random sampling.
#'  
#' @details 
#'  Calculates exact p-values for permutation tests when permutations are 
#'  randomly drawn with replacement using the \code{\link[statmod]{permp}} 
#'  function in the \code{\link{statmod}} package.
#'  
#' @references 
#'   \enumerate{
#'     \item{
#'       Phipson, B. & Smyth, G. K. \emph{Permutation P-values should never be 
#'       zero: calculating exact P-values when permutations are randomly drawn.}
#'       Stat. Appl. Genet. Mol. Biol. \strong{9}, Article39 (2010). 
#'     }
#'   }
#'  
#' @param permuted vector of values making up the empirical distribution.
#' @param observed the observed value of the test statistic.
#' @param subsetSize the size of the network subset the null distribution is 
#'  drawn for.
#' @param totalSize the size of the whole network
#' @param order logical; does the order of nodes in the permutation affect the 
#'  value of the test statistic?
#' @param alternative a character string specifying the alternative hypothesis, 
#'  must be one of "greater" (default), "less", or "two.sided". 
#'  You can specify just the initial letter.
#'  
#' @aliases permutation permuted
#' @name permutation
#' @export
perm.test <- function(
  permuted, observed, subsetSize, totalSize, order=TRUE, alternative="greater"
) {
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  
  if (is.na(observed))
    return(NA)
  
  if (order) {
    total.nperm = prod(totalSize:(totalSize - subsetSize + 1))
  } else {
    total.nperm = choose(totalSize, subsetSize)
  }
  permuted <- sort(permuted)
  nPerm <- length(permuted)
  
  less.extreme <- length(permuted[permuted <= observed])
  more.extreme <- length(permuted[permuted >= observed])
  lower.pval <- permp(less.extreme, nPerm, total.nperm=total.nperm)
  upper.pval <- permp(more.extreme, nPerm, total.nperm=total.nperm)
  
  if (altMatch == 1L) {
    return(min(lower.pval, upper.pval)*2)
  } else if (altMatch == 2L) {
    return(lower.pval) 
  } else if (altMatch == 3L) {
    return(upper.pval)
  } 
}

#' Exact permutation p-values wrapper
#' 
#' Wrapper for \code{\link[statmod]{permp}} from the 
#' \code{\link[statmod]{statmod}} library, which can crash if FORTRAN 
#' libraries are not properly linked.
#' 
#' @details
#' In the case \code{\link[statmod]{permp}} fails, the wrapper will fall back 
#' to a slightly more conservative biased estimator: (1+x)/(1+nPerm).
#' 
#' @param x number of permutations that yielded test statistics at least as 
#'  extreme as the observed data. May be a vector or an array of values. 
#' @param nperm total number of permutations performed.
#' @param ... other arguments to pass to\code{\link[statmod]{permp}}.
#' @return 
#'  vector or array of p-values, of same dimensions as \code{x}.
#' @importFrom statmod permp
permp <- function(x, nperm, ...) {
  tryCatch({
    return(statmod::permp(x, nperm, ...))
  }, error=function(e) {
    warning(
      "Error from statmod::permp:", e$message, 
      "\nUsing conservative biased estimator (1+x)/(1+nPerm) instead."
    )
    return(
      (x + 1)/(nperm + 1)  
    )
  })
}

#' @description
#'  \code{requiredPerms}:  how many permutations do I need to be able to detect
#'  significance at a given threshold \code{alpha}?
#' 
#' @param alpha desired significance threshold.
#' @return The minimum number of permutations required to detect any significant
#'  associations at the provided \code{alpha}. The minimum p-value will always
#'  be smaller than \code{alpha}.
#' @rdname permutation
#' @export
requiredPerms <- function(alpha, alternative="greater") {
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  
  if (altMatch == 1) {
    1/alpha*2
  } else {
    1/alpha
  }
}
