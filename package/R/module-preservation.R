# This is the core functionality of the package.

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
#' @param ignoreCol A vector of \code{moduleLabels} colors to ignore when
#'  calculating preservation.
#' @param verbose Logical; If \code{TRUE}, outputs message to file.
#' @param logFile File to send output to, if \code{verbose} is \code{TRUE}.
#' @param savePermuted
#' @param permutedFile
#' @return a
#' 
#' @export
modulePreservation <- function(expressionSets=NULL, adjacencySets=NULL, 
                               referenceSets, testSets, moduleLabels, 
                               ignoreCol = NULL, nPermutations=200, 
                               verbose=TRUE, logFile='log.txt', save=TRUE,
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
      # Restrict to overlapping genes between reference and test.
      if (!is.null(adjacencySets)) {
        refGenes <- rownames(adjacencySets[[ref]])
        testGenes <- rownames(adjacencySets[[test]])
      } else {
        refGenes <- rownames(expressionSets[[ref]])
        testGenes <- rownames(expressionSets[[test]])
      }
      overlapGenes <- intersect(refGenes, testGenes)
      if(length(overlapGenes) == 0) {
        warning("No overlap of genes between reference network", ref, 
                "and test network", test, ". Skipping to next pair.")
        return(NULL)
      } else if (verbose) {
        cat("Calculating preservation on overlap of", length(overlapGenes), "/",
            length(refGenes), "genes between network pair.\n")
      }
      
      # Get Information about the modules
      moduleSizes <- table(moduleLabels[[ref]])
      overlapModSizes <- table(moduleLabels[[ref]][overlapGenes])
      modOverlap <- overlapModSizes/moduleSizes
      
      moduleNames <- unique(moduleLabels[[ref]][overlapGenes])
      moduleNames <- moduleNames[!(moduleNames %in% ignoreCol)]
      
      # Create dictionary of module names -> nodes.
      moduleIndices <- sapply(moduleNames, function(name) {
        which(moduleLabels[[ref]][overlapGenes] == name)
      })
      
      # Calculate observed preservation statistics for the modules.
      cat("Calculating observed statistics...\n")
      observed <- foreach(module=moduleNames, .combine=rbind) %do% {
        calculatePreservation(expressionSets[[ref]], adjacencySets[[ref]],
                              expressionSets[[test]], adjacencySets[[test]],
                              moduleIndices[[module]])
      }
      rownames(observed) <- moduleNames
      
      # Calculate the null distribution for each of the preservation statistics.
      permuted <- foreach(i=1:nPermutations, .combine=.bind3) %dopar% {
        cat("Permutation ", i, "\n")
        thisPerm <- foreach(module=moduleNames, .combine=rbind) %do% {
          # Randomise module assignments for each run
          permuted <- sample(1:length(overlapGenes), size=moduleSizes[module])
          
          calculatePreservation(expressionSets[[ref]], adjacencySets[[ref]],
                                expressionSets[[test]], adjacencySets[[test]],
                                permuted)
        }
        rownames(thisPerm) <- moduleNames
        return(thisPerm)
      }
      
      # Calculate the p-value for the observed statistic based on the NULL 
      # distribution
      cat("Calculating P-values...")
      p.values <- foreach(j=1:ncol(observed), .combine=cbind) %:% 
                    foreach(i=1:nrow(observed), .combine=c) %do% {
        pperm(permuted[i,j,], observed[i,j], lower.tail=FALSE)
      }
      dimnames(p.values) <- dimnames(observed)

      return(list(observed=observed, permuted=permuted, p.value=p.values,
                  overlap=modOverlap[moduleNames], 
                  overlapSize=overlapModSizes[moduleNames]))
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
#' @param moduleNodes Nodes (genes) in the reference network for the module 
#'    whose preservation is being measured.
#' @return A vector of preservation scores for the module.
#'   
calculatePreservation <- function(refExpr, refAdj, testExpr, testAdj,
                                  moduleNodes) {

  # if provided, must be for both reference and test networks.
  stopifnot(!xor(is.null(refExpr), is.null(testExpr)))
  stopifnot(!xor(is.null(refAdj), is.null(testAdj)))
  
  exprPres <- NULL
  adjPres <- NULL
  if (!is.null(refExpr)) {
    
  }  
  if (!is.null(refAdj)) {
    adjNames <- c("meanAdj", "corKIM", "corKALL")
    adjPres <- c(
      meanAdj(testAdj, moduleNodes),
      cor(kIM(refAdj, moduleNodes, FALSE), kIM(testAdj, moduleNodes, FALSE), method="spearman"),
      cor(kIM(refAdj, moduleNodes, TRUE), kIM(testAdj, moduleNodes, TRUE), method="spearman")
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
