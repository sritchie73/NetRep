# This is the core functionality of the package.

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
#' @param moduleLabels
#' @param verbose
#' @param savePermuted
#' @param permutedFile
#' @return a
#' 
#' @import foreach
#' @import abind
#' @export
modulePreservation <- function(expressionSets=NULL, adjacencySets=NULL, 
                               referenceSets, testSets, moduleLabels,
                               nPermutations=200, verbose=TRUE, save=TRUE,
                               permutedFile="permutedStatistics.rda") {
  # TODO: Check input arguments for errors and consistency
  if (is.null(expressionSets) && is.null(adjacencySets)) {
    stop("Nothing to calculate!")
  } 
  
  nNets <- length(adjacencySets)
  
  # Iterate over networks pairwise, and calculate preservation between reference
  # and test pairs.
  preservation <- foreach (ref=1:nNets) %:% foreach(test=1:nNets) %do% {
    if ((ref %in% referenceSets) && (test %in% testSets)) {
      cat("Calculating preservation for network pair ", ref, ", ", test, "\n")
      # Get Information about the modules
      # TODO: restrict to overlapping genes
      moduleSizes <- table(moduleLabels[[ref]])
      moduleNames <- unique(moduleLabels[[ref]])
      
      # TODO: Handle where no adjnet provided
      nRefNodes <- ncol(adjacencySets[[ref]])
      nTestNodes <- ncol(adjacencySets[[test]])
      
      # Create dictionary of module names -> nodes.
      moduleIndices <- sapply(moduleNames, function(name) {
        which(moduleLabels[[ref]] == name)
      })
      
      # Calculate observed preservation statistics for the modules.
      # TODO: Handle modIndexing of testNetwork properly
      cat("Calculating observed statistics...\n")
      observed <- foreach(module=moduleNames, .combine=rbind) %do% {
        calculatePreservation(expressionSets[[ref]], adjacencySets[[ref]],
                              expressionSets[[test]], adjacencySets[[test]],
                              moduleIndices[[module]], moduleIndices[[module]])
      }
      rownames(observed) <- moduleNames
      
      # Calculate the null distribution for each of the preservation statistics.
      permuted <- foreach(i=1:nPermutations, .combine=.bind3) %dopar% {
        cat("Permutation ", i, "\n")
        thisPerm <- foreach(module=moduleNames, .combine=rbind) %do% {
          # Generate permutation indices for this module
          permutedIndices <- sample(1:nRefNodes, size=moduleSizes[module])
          
          calculatePreservation(expressionSets[[ref]], adjacencySets[[ref]],
                                expressionSets[[test]], adjacencySets[[test]],
                                moduleIndices[[module]], permutedIndices)
        }
        rownames(thisPerm) <- moduleNames
        return(thisPerm)
      }
      
      # Calculate the p-value for the observed statistic based on the NULL 
      # distribution
      p.values <- foreach(j=1:ncol(observed), .combine=cbind) %:% 
                    foreach(i=1:nrow(observed), .combine=c) %do% {
        pperm(permuted[i,j,], observed[i,j], lower.tail=FALSE)
      }
      dimnames(p.values) <- dimnames(observed)

      return(list(observed=observed, permuted=permuted, p.value=p.values))
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
#' @param testAdj Adjacency matrix for the test network.
#' @param refModuleNodes Nodes (genes) in the reference network for the module 
#'    whose preservation is being measured.
#' @param testModuleNodes Nodes (genes) in the test network for the module whose 
#'    preservation is being measured.
#' @return A vector of preservation scores for the module.
#'   
calculatePreservation <- function(refExpr, refAdj, testExpr, testAdj,
                                  refModuleNodes, testModuleNodes) {
  # Basic Error Checking
  stopifnot(length(refModuleNodes) == length(testModuleNodes))
  # if provided, must be for both reference and test networks.
  stopifnot(!xor(is.null(refExpr), is.null(testExpr)))
  stopifnot(!xor(is.null(refAdj), is.null(testAdj)))
  
  exprPres <- NULL
  adjPres <- NULL
  if (!is.null(refExpr)) {
    
  }  
  if (!is.null(refAdj)) {
    adjNames <- c("meanAdj", "meanAdj2")
    adjPres <- c(
      meanAdj(testAdj, refModuleNodes, FALSE),
      meanAdj(testAdj, refModuleNodes, FALSE)
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
