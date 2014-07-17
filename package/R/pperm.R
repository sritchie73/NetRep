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
#'  to estimate the p-value of \code{q} by setting \code{tailApprox} to 
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
#'     \item{
#'       Knijnenburg, T. A., Wessels, L. F. A., Reinders, M. J. T. & Shmulevich, 
#'       I. \emph{Fewer permutations, more accurate P-values}. Bioinformatics 
#'       \strong{25}, i161-8 (2009). 
#'     }
#'   }
#'  
#' @param permuted vector of values making up the empirical distribution.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is 
#'  taken to be the number required. 
#' @param log.p logicial; if TRUE, probabilities p are given as log(p)
#' @param tailApprox logical; if \code{TRUE}, use the tail approximation 
#'  algorithm to estimate extreme p-values (see details).
#' @param lower.tail logical; if TRUE (default), probabilities are 
#'    \eqn{P[x \le x]} otherwise \eqn{P[X > x]}.   
#'  
#' @examples
#'  # A contrived example. For large n these results will be the same as qnorm, 
#'  # rnorm, pnorm.
#'  normData <- rnorm(n=10000)
#'  pperm(normData, -1.644854) # should be approximately 0.05
#'  qperm(normData, 0.95) # should be approximately 1.644854
#'  rperm(normData, 100) # should be similar to rnorm(100)
#'    
#'  
#' @aliases permutation permuted
#' @name permutation
NULL
  
#' @rdname permutation
#' @export
#' @import gPdtest
#' @importFrom fExtremes pgpd
pperm <- function(permuted, q, tailApprox=FALSE, lower.tail=TRUE, log.p=FALSE) {
  permuted <- sort(permuted)
  if (lower.tail) {
    more.extreme <- length(permuted[permuted < q])
  } else {
    more.extreme <- length(permuted[permuted > q])
  }
  if (tailApprox & more.extreme < 10) {
    # Starting with the 250 most extreme observations (exceedances) , decrease 
    # by 10, until we have a good fit to a GPD (p > 0.05), see (2).
    topn <- 250
    while (topn > 0) {
      test <- gpd.test(tail(permuted, topn))
      ifelse(test$boot.test$p.value > 0.05, break, topn <- topn - 10)
    }
    if (topn == 0) {
      warning("unable to fit a generalized pareto distribution to the permuted",
              " data!")
      p.value <- pperm(permuted, q, FALSE, lower.tail, log.p)
    } else {
      fit <- gpd.fit(tail(permuted, topn), method="amle")
      p.value <- pgpd(q, xi=fit[[1]], beta=fit[[2]], lower.tail=lower.tail)
      p.value <- as.numeric(p.value) # get rid of extra attributes
      if (log.p) {
        p.value <- log(p.value)
      }
    }
  } else {
    # If we're using tail approximation, no need for the pseudocount (2)
    pCount <- ifelse(tailApprox, 0, 1)
    nPerm <- length(permuted)
    if (log.p) {
      p.value <- log(more.extreme + pCount) - log(nPerm + pCount)
    } else {
      p.value <- (more.extreme + pCount) / (nPerm + pCount)
    }
  }
  return(p.value)
}

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
  # Index/ices of the closest p-values to the given one
  p.index <- which(is.equal(p.diff, min(p.diff)))
  if (length(p.index) > 1) {
    warning("Specified quantile ", p, " is halfway between two points on the",
            " permuted distribution!")
  }
  return(sort(permuted)[p.index])
}


#' @usage rperm(permuted, n)
#' @rdname permutation
#' @export
rperm <- function(permuted, n) {
  sample(permuted, size=n, replace=TRUE)
}

#' @description
#'  \code{requiredPower}: If the tail approximation is not being used, how many
#'  permutations do I need to be able to detect significance at a given
#'  threshold \code{alpha}?
#' 
#' @param alpha desired significance threshold.
#' @return The minimum number of permutations required to detect any significant
#'  associations at the provided \code{alpha}.
#' @rdname permutation
#' @export
requiredPower <- function(alpha) {
  1/(alpha) - 1
}

#' Floating Point Comparison
#'
#' Tests elements in a vector for equality to a specified floating point value.
#' 
#' @param vector a numeric vector of doubles
#' @param value a double
#' @return logical; \code{TRUE} where an element is equal to \code{value}, 
#'   \code{NA} where a comparison is made with an \code{NA}, and \code{FALSE}
#'   otherwise.
is.equal <- function(vector, value) {
  sapply(vector, function(element) {
    isTRUE(all.equal(element, value))
  })
}


