# This is the core functionality of the package.

library(bigmemory)
library(biganalytics)
library(foreach)
library(abind)

#' Calculates pairwise between reference and test networks 
#' 
#' Description goes here
#' 
#' @param expressionSets
#' @param adjacencySets
#' @param referenceSets 
#' @param testSets
#' @param nPermutations
#' @param refModLabels
#' @param verbose
#' @param savePermuted
#' @param permutedFile
#' @return a
#' @export
modulePreservation <- function(expressionSets=NULL, adjacencySets=NULL, 
                               referenceSets, testSets, refModLabels,
                               nPermutations=200, verbose=TRUE, save=TRUE,
                               permutedFile="permutedStatistics.rda") {
  # TODO: Check input arguments for errors and consistency
  if (is.null(expressionSets) && is.null(adjacencySets)) {
    stop("Nothing to calculate!")
  } 
  
  nNets <- length(adjacencySets)
  
  # Change diagonals of Adjacency to NA for preservation statistics:
  if (!is.null(adjacencySets)) {
    lapply(adjacencySets, function(adj) { 
      diag(adj) <- NA 
    })
  }
  
  # Convert to big.matrix descriptors so that the networks only have to be 
  # loaded in memory once across all threads.
  # expression/adjacencySets hold the actual big.matrix
  # expression/adjacencyDesc hold the pointers for the threads to use.
  # If i don't keep the Sets, the matrices will be garbage collected.
  expressionDesc <- NULL
  adjacencyDesc <- NULL
  if (!is.null(expressionSets)) {
    expressionSets <- lapply(expressionSets, .toBigMatrix)
    expressionDesc <- lapply(expressionSets, describe)
  }
  if (!is.null(adjacencySets)) {
    adjacencySets <- lapply(adjacencySets, .toBigMatrix)
    adjacencyDesc <- lapply(adjacencySets, describe)
  }
  
  # Iterate over networks pairwise, and calculate preservation between reference
  # and test pairs.
  preservation <- foreach (ref=1:nNets) %:% foreach(test=1:nNets) %do% {
    if ((ref %in% referenceSets) && (test %in% testSets)) {
      # Get Information about the modules
      # TODO: restrict to overlapping genes
      modSizes <- table(refModLabels[[ref]])
      modNames <- unique(refModLabels[[ref]])
      modIndexes <- sapply(modNames, function(mod) { # hash by module name
        which(modNames == mod)
      })
      
      # TODO: Get the observed statistics for each of the modules
      foreach(j=1:(length(modNames) + 1), .combine=rbind) %do% {
        
      }
      
      
      # TODO: Generate the null distribution for each statistic
      permuted <- foreach(i=1:nPermutations, .combine=.bind3) %dopar% {
        require(bigmemory)
        
      }
      
      # TODO: Calculate p-value for each statistic
      
    } else {
      return(NULL)
    }
  }
  return(preservation)
}

#' Calculates the preservation of a module between two networks using a wide
#' range of statistics.
#' 
#' @details This is the actual guts of the calculation. This calculates the 
#' preservation of each module based on the statistics outlined by Langfelder
#' et al. (See Reference below) and implemented in their R package WGCNA.
#' For each module the following are calculated depending on whether the 
#' expression data, adjacency data, or both are provided:
#' \itemize{
#'  \item{Expression Data}{
#'    \itemize{
#'      \item{TODO:}{None Implemented yet}
#'    }
#'  }
#'  \item{Network Adjacencies}{
#'    \itemize{
#'      \item{\code{\link{meanAdj}}:}{
#'        The average edge weight of the module in the test network.
#'      }
#'    }
#'  }
#' }
#'
#' @references Langfelder, Peter; Luo, Rui; Oldham, Michael C.; and Horvath, 
#'    Steve. Is My Network Module Preserved and Reproducible?. PLoS 
#'    Computational Biology, 2011.
#' 
#' @param refExpr Expression matrix for the reference network.
#' @param refAdj Adjacency matrix for the reference network.
#' @param testExpr Expression matrix for the test network.
#' @param refModNodes Nodes (genes) in the reference network for the module 
#'    whose preservation is being measured.
#' @param testModNodes Nodes (genes) in the test network for the module whose 
#'    preservation is being measured.
#' @return A vector of preservation scores for the module.
#'   
calculatePreservation <- function(refExpr, refAdj, testExpr, testAdj,
                                   refModNodes, testModNodes) {
  # Basic Error Checking
  stopifnot(length(refNodes) == length(testNodes))
  # if provided, must be for both reference and test networks.
  stopifnot(!xor(is.null(refExpr), is.null(testExpr)))
  stopifnot(!xor(is.null(refAdj), is.null(testAdj)))
  
  exprPres <- NULL
  adjPres <- NULL
  if (!is.null(refExpr)) {
 
  }
  if (!is.null(refAdj)) {
    adjNames <- c("meanAdj")
    adjPres <- c(
      meanAdj(testAdj, testModNodes)
    )
    names(adjPres) <- adjNames
  }
  preservation <- c(exprPres, adjPres)
  return(preservation)
}

#####################
#
# Helper functions
#
#####################

# Bind two dimensional objects along a third dimension
.bind3 <- function(...) { abind(..., along=3) }

# Custom converter
.toBigMatrix <- function(mat) { 
  if(!is.big.matrix(mat)) {
    return(as.big.matrix(mat))
  }
}
