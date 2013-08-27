#' @title Distribution from Permutations
#' 
#' @description Distribution function, quantile function, and random generation
#'  for an empirical distribution provided by the user in the \code{permuted}
#'  argument. These distributions are typically drawn through permutation 
#'  testing of some statistic.
#'  
#' @rdname permutation
#' @references Phipson, Belinda; and Smyth, Gordon. Permutation p-values should
#'  never be zero: calculating exact p-values when permutations are randomly
#'  drawn. Statistical Applications in Genetics and Molecular Biology (9), 2010.
#'  
#' @param permuted vector of values making up the empirical distribution.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is 
#'  taken to be the number required. 
#' @param log.p logicial; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are 
#'    \eqn{P[x \le x]} otherwise \eqn{P[X > x]}.   
#'  
#' @examples
#' # A contrived example. For large n these results will be the same as qnorm, 
#' # rnorm, pnorm.
#' normData <- rnorm(n=10000)
#' pperm(normData, -1.644854) # should be approximately 0.05
#' qperm(normData, 0.95) # should be approximately 1.644854
#' rperm(normData, 100) # should be similar to rnorm(100)
#'    
#'  
#' @aliases pperm, qperm, rperm, permutation, permuted
#'  


#' @details P-values are calculated by \code{pperm} using proportions on the 
#'  provided distribution (\code{permuted}). An observation, \code{q}'s, p-value 
#'  is the proportion of data from \code{permuted} whom are more extreme than 
#'  \code{q}. This is accurate strictly when:
#'  \enumerate{
#'    \item Each permutation is independent.
#'    \item Permutations are sampled without replacement.
#'  }
#'  For our purposes this calculation is a good enough approximation, but when 
#'  more accuracy is required see the package \code{rperm}.
#'  
#' @usage pperm(permuted, q, lower.tail = TRUE, log.p = FALSE) 
#' @rdname permutation
#' @export
pperm <- function(permuted, q, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail) {
    more.extreme <- length(permuted[permuted < q])
  } else {
    more.extreme <- length(permuted[permuted > q])
  }

  nPermutations <- length(permuted)
  # Addition of 1 to both numerator and denominator makes it impossible to get
  # a p-value of 0. See the @reference reading.
  if (log.p) {
    p.value <- log(more.extreme + 1) - log(nPermutations + 1)
  } else {
    p.value <- {more.extreme + 1} / {nPermutations + 1}
  }
  return(p.value)
}

#' @details It is rarely possible to give an accurate value for \code{qperm} 
#'  since there is no theoretical distribution to draw from. Instead the closest
#'  possible observation from the data is returned.
#'  NOTE: It is possible for the quantile \code{p} to fall exactly halfway 
#'  between two observations from the \code{permuted} distribution. In this case
#'  both data points are returned and a warning is thrown.
#'  It is up to the user to choose which observation to take, the conservative
#'  approach would be to choose the first observation if \code{p} < 0.5, or the
#'  second if \code{p} > 0.5. For large \code{permuted} distributions the 
#'  interpretation will not change much.
#' 
#' @usage qperm(permuted, p, log.p = FALSE) 
#' @rdname permutation
#' @export
qperm <- function(permuted, p, log.p = FALSE) {
  if (((p < 0) || (p > 1)) && (!log.p)) {
    stop("Error: `p` must be between 0 and 1.")
  } else if ((p > 0) && (log.p)) {
    stop("Error: `log(p)` must be between -Inf and 0.")
  }
  
  # Generate all possible quantiles and pick the closest one
  nPermutations <- length(permuted)
  if (log.p) {
    p.possible <- log(0:(nPermutations - 1) + 1) - log(nPermutations + 1)
  } else {
    p.possible <- {0:(nPermutations - 1) + 1} / {nPermutations + 1} 
  }
  p.diff <- abs(p.possible - p)
  p.closest <- p.possible[which.equal(p.diff, min(p.diff))]
  if (length(p.closest) > 1) {
    warning("Specified quantile ", p, " is halfway between two points on the",
            " permuted distribution!")
  }
  p.index <- p.closest*(nPermutations + 1)
  return(sort(permuted)[p.index])
}


#' @usage rperm(permuted, n)
#' @rdname permutation
#' @export
rperm <- function(permuted, n) {
  sample(permuted, size=n, replace=TRUE)
}

#' Finds elements in a vector which are equal to a specified floating point 
#' value
which.equal <- function(vector, value) {
  sapply(vector, function(element) {
    isTRUE(all.equal(element, value))
  })
}
