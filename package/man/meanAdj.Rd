\name{meanAdj}
\alias{meanAdj}
\title{Calculate the meanAdj statistic of a module in a test network}
\usage{
  meanAdj(testNet, testIndices)
}
\arguments{
  \item{testNet}{the adjacency matrix for the test
  network.}

  \item{testIndices}{the indices on the nodes in the
  module.}
}
\value{
  The mean edge weight of the module in the test network.
}
\description{
  Calculates the mean edge weight of the module in the test
  network. A high preservation score indicates preservation
  of the module's density.
}
\note{
  Make sure to mark with NA edges that should not be
  included in the calculation. This is done to minimise
  computation. For example, gene coexpression adjacency
  matrices from WGCNA should have their diagonal replaced
  with NA. For sparse networks the user will need to
  determine where NA or 0 is most appropriate.
}
\references{
  Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and
  Horvath, Steve. Is My Network Module Preserved and
  Reproducible?. PLoS Computational Biology, 2011.
}

