\name{meanAdj}
\alias{meanAdj}
\title{Calculate the meanAdj statistic of a module in a test network}
\usage{
  meanAdj(adjacency, moduleIndices, includeDiagonal)
}
\arguments{
  \item{adjacency}{Adjacency matrix to calculate meanAdj of
  module in.}

  \item{moduleIndices}{Where the nodes in the module are
  located in the adjacency network.}

  \item{includeDiagonal}{Logical, if FALSE, diagonal values
  (edges between nodes and themselves) are not counted
  towards the mean edge weight.}
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

