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
#'  must be one of "greater" (default), "less", or "two.sider". 
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
  
  if (order) {
    total.nperm = prod(totalSize:(totalSize - subsetSize + 1))
  } else {
    total.nperm = choose(totalSize, subsetSize)
  }
  permuted <- sort(permuted)
  nPerm <- length(permuted)
  
  less.extreme <- length(permuted[permuted < observed])
  more.extreme <- length(permuted[permuted > observed])
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
#'  \code{requiredPower}:  how many permutations do I need to be able to detect
#'  significance at a given threshold \code{alpha}?
#' 
#' @param alpha desired significance threshold.
#' @return The minimum number of permutations required to detect any significant
#'  associations at the provided \code{alpha}. The minimum p-value will always
#'  be smaller than \code{alpha}.
#' @rdname permutation
#' @export
requiredPower <- function(alpha, alternative="greater") {
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


#' Visual Sanity Check of pperm
#' 
#' Generate a qqplot of the observed p-values from \code{\link{pperm}}, vs the
#' p-values we would expect from the theoretical distribution the 
#' \code{permuted} observations were drawn from (given in \code{tfunc}).
#' 
#' @note
#'  For internal testing purposes.
#' 
#' @param permuted vector of values making up the empirical distribution.
#' @param q vector of quantiles to test.
#' @param tfunc a probability function to use to get the expected p-values. This
#'  should correspond to the theortical distribution \code{permuted} was drawn
#'  from.
#' @param tail.approx logical; if \code{TRUE}, use the tail approximation 
#'  algorithm to estimate extreme p-values (see details).
#' @param lowerTail logical; if TRUE (default), probabilities are 
#'  \eqn{P[x \le x]} otherwise \eqn{P[X > x]}.
#' @param ... other arguments to pass to \code{\link[graphics]{plot}}. Note that
#'   providing \code{xlab}, \code{ylab}, \code{xlim}, or \code{ylim} will result
#'   in an error. 
#' @importFrom graphics plot
#' @export
qqperm <- function(permuted, q, tfunc, tail.approx=TRUE, lowerTail=FALSE, ...) {
  permuted <- sort(permuted)
  expected <- do.call(tfunc, list(q=q, lowerTail=lowerTail))
  observed <- sapply(q, function(x) {
    pperm(permuted, x, lowerTail=lowerTail, tail.approx=tail.approx)
  })
  min <- min(c(-log10(expected), -log10(observed)))
  max <- max(c(-log10(expected), -log10(observed)))
  plot(-log10(expected), -log10(observed), xlab="expected", ylab="observed",
       xlim=c(min, max), ylim=c(min, max), ...)
  abline(0, 1)
  # Color p-values that have been estimated using the tail approximation method
  if (tail.approx) {
    if (lowerTail) {
      tai <- which(q < permuted[11])
    } else {
      tai <- which(q > tail(permuted, 11)[1])
    }
    points(-log10(expected[tai]), -log10(observed[tai]), col="red")
  } 
}
