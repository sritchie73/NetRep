\name{doParPermute}
\alias{doParPermute}
\title{Run permutations on the provided function}
\usage{
  doParPermute(FUN, permuteVector, nPerm = 1000,
    verbose = TRUE, logFile = "log.txt", ...)
}
\arguments{
  \item{FUN}{the function to pass the \code{permuteVector}
  to.}

  \item{permuteVector}{vector to shuffle. Names of the
  vector will remain constant for each permutation.}

  \item{nPerm}{Number of permutations to run}

  \item{verbose}{Logical; If TRUE, writes the current
  progress to file at every iteration.}

  \item{logFile}{The file to output progress to if
  \code{verbose} is \code{TRUE}.}

  \item{...}{optional arguments to \code{\link{foreach}}}
}
\description{
  Does \code{nPerm} calculations on the provided expression
  in parallel. The \code{permuteVector} is shuffled, and
  passed to the provided function.
}

