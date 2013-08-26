\name{meanAdj}
\alias{meanAdj}
\title{Calculate the meanAdj statistic of a module in a test network}
\usage{
  meanAdj(testNet, testIndices, diagonal = NA)
}
\arguments{
  \item{testNet}{the adjacency matrix for the test
  network.}

  \item{testIndices}{the indices on the nodes in the
  module.}

  \item{diagonal}{How to handle diagonals in the adjacency
  matrix: \itemize{ \item{NA }{An edge from a node to
  itself is allowed and counted in the mean.} \item{numeric
  }{The value the all diagonal entries of the adjacency
  matrix have, to be removed in the calculation of the
  mean.} }}
}
\value{
  The mean edge weight of the module in the test network.
}
\description{
  Calculates the mean edge weight of the module in the test
  network. A high preservation score indicates preservation
  of the module's density.
}
\references{
  Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and
  Horvath, Steve. Is My Network Module Preserved and
  Reproducible?. PLoS Computational Biology, 2011.
}

