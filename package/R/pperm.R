#' @title Distribution from Permutations
#' 
#' @description Distribution function, quantile function, and random generation
#'  for an empirical distribution provided by the user in the \code{permuted}
#'  argument. These distributions are typically drawn through permutation 
#'  testing of some statistic.
#'  
#' @details 
#'  P-values are calculated by \code{pperm} using proportions on the provided 
#'  distribution (\emph{permuted}). An observation, \code{q}'s, p-value 
#'  is the proportion of \emph{permuted} whom are more extreme than \code{q}. 
#'  This is accurate strictly when:
#'  \enumerate{
#'    \item Each permutation is independent.
#'    \item Permutations are sampled without replacement.
#'  }
#'  While these assumptions may not necessarily always be true, for our purposes
#'  they are a reasonable approximation. If \code{q} is more extreme than the 
#'  \code{permuted}, the tail approximation method from \emph{(2)} can be used
#'  to estimate the p-value of \code{q} by setting \code{tail.approx} to 
#'  \code{TRUE}.
#'  
#'  It is rarely possible to give an accurate value for \code{qperm} since there
#'  is no theoretical distribution to draw from. Instead the closest possible
#'  observation from the data is returned.
#'  
#' @note
#'  It is possible for the quantile \code{p} to fall exactly halfway 
#'  between two observations from the \emph{permuted} distribution. In this case
#'  both data points are returned and a warning is generated. It is up to the 
#'  user to choose which observation to take: the conservative approach is to 
#'  choose the first observation if \code{p} < 0.5, or the second observation if
#'  \code{p} > 0.5.
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
#' @param subsetSize the size of the network subset the null distribution is 
#'  drawn for.
#' @param q vector of quantiles.
#' @param lowerTail logical; if TRUE (default), probabilities are 
#'    \eqn{P[x \le x]} otherwise \eqn{P[X > x]}.   
#'  
#' @examples
#'  # A contrived example. For large n these results will be the same as qnorm, 
#'  # rnorm, pnorm.
#'  normData <- rnorm(n=10000)
#'  pperm(normData, 1.644854, 1) # should be approximately 0.05
#'  
#' @aliases permutation permuted
# @param tail.approx logical; if \code{TRUE}, use the tail approximation 
#  algorithm to estimate extreme p-values (see details).
#' @name permutation
#' @importFrom statmod permp
#' @export
pperm <- function(permuted, q, moduleSize, lowerTail=FALSE) {
  permuted <- sort(permuted)
  if (lowerTail) {
    more.extreme <- length(permuted[permuted < q])
  } else {
    more.extreme <- length(permuted[permuted > q])
  }
  nPerm <- length(permuted)
  p.value <- permp(more.extreme, nPerm, total.nperm=choose(nPerm, moduleSize))
  return(p.value)
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
requiredPower <- function(alpha) {
  1/alpha
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
