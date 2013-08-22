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
      # TODO: Get the genes in the reference modules which are present in the 
      #   Test network.
      
      # TODO: Get the observed statistics for each of the modules
      
      
      # TODO: Generate the null distribution for each statistic
      permuted <- foreach(1:nPermutations, .combine=.bind3) %dopar% {
        require(bigmemory)
        
      }
    } else {
      return(NULL)
    }
  }
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
