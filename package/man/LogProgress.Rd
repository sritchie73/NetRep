\name{LogProgress}
\alias{LogProgress}
\title{Log the progress of a looping function to file.}
\usage{
  LogProgress(this, total, file = "log.txt")
}
\arguments{
  \item{this}{the current iteration number.}

  \item{total}{the total number of iterations}

  \item{file}{The log file to output to.}
}
\description{
  Outputs the progress of a looping function to a file.
  This can then be monitored using the unix command
  \code{tail -f log.txt}. Useful for parallel loops, where
  \code{\link{PrintProgress}} won't be able to output to
  the R terminal
}

