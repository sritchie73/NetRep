\name{pperm}
\alias{permutation,}
\alias{permuted}
\alias{pperm}
\alias{pperm,}
\alias{qperm}
\alias{qperm,}
\alias{rperm}
\alias{rperm,}
\title{Distribution from Permutations}
\usage{
  pperm(permuted, q, lower.tail = TRUE, log.p = FALSE)

  qperm(permuted, p, log.p = FALSE)

  rperm(permuted, n)
}
\arguments{
  \item{permuted}{vector of values making up the empirical
  distribution.}

  \item{q}{vector of quantiles.}

  \item{p}{vector of probabilities.}

  \item{n}{number of observations. If \code{length(n) > 1},
  the length is taken to be the number required.}

  \item{log.p}{logicial; if TRUE, probabilities p are given
  as log(p)}

  \item{lower.tail}{logical; if TRUE (default),
  probabilities are \eqn{P[x \le x]} otherwise \eqn{P[X >
  x]}.}
}
\description{
  Distribution function, quantile function, and random
  generation for an empirical distribution provided by the
  user in the \code{permuted} argument. These distributions
  are typically drawn through permutation testing of some
  statistic.
}
\details{
  P-values are calculated by \code{pperm} using proportions
  on the provided distribution (\code{permuted}). An
  observation, \code{q}'s, p-value is the proportion of
  data from \code{permuted} whom are more extreme than
  \code{q}. This is accurate strictly when: \enumerate{
  \item Each permutation is independent.  \item
  Permutations are sampled without replacement.  } For our
  purposes this calculation is a good enough approximation,
  but when more accuracy is required see the package
  \code{rperm}.

  It is rarely possible to give an accurate value for
  \code{qperm} since there is no theoretical distribution
  to draw from. Instead the closest possible observation
  from the data is returned.  NOTE: It is possible for the
  quantile \code{p} to fall exactly halfway between two
  observations from the \code{permuted} distribution. In
  this case both data points are returned and a warning is
  thrown.
}
\examples{
# A contrived example. For large n these results will be the same as qnorm,
# rnorm, pnorm.
normData <- rnorm(n=10000)
pperm(normData, -1.644854) # should be approximately 0.05
qperm(normData, 0.95) # should be approximately 1.644854
rperm(normData, 100) # should be similar to normData and rnorm(100)
}
\references{
  Phipson, Belinda; and Smyth, Gordon. Permutation p-values
  should never be zero: calculating exact p-values when
  permutations are randomly drawn. Statistical Applications
  in Genetics and Molecular Biology (9), 2010.
}

