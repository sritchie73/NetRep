#' Calculate module preservation statistics.
#' 
#' Calculates the modules preservation statistics using the pre-calculated 
#' module properties along with the correlation structure matrices.
#' 
#' @param discProps module properties in the \emph{discovery} dataset calculated 
#'  using \code{\link{netProps}} and \code{\link{dataProps}}.
#' @param testProps module properties in the \emph{test} dataset calculated 
#'  using \code{\link{netProps}} and \code{\link{dataProps}}.
#' @param discCor 'bigMatrix' containing the correlation structure matrix 
#'  for the \emph{discovery} dataset.
#' @param discIndices module indices in the discovery dataset.
#' @param testCor 'bigMatrix' containing the correlation structure matrix 
#'  for the \emph{test} dataset.
#' @param testIndices module indices in the test dataset.
#'  
#' @section Design decision:
#'  Most of the module presevation statistics are calculated as the correlation 
#'  of some network property between the \emph{discovery} and \emph{test} 
#'  datasets. The permutation procedure randomly samples from the \emph{test} 
#'  dataset, so it therefore makes sense to calculate the network properties in 
#'  the \emph{discovery} dataset once only. The code is thus split into 
#'  calculations of the network properties, and then subsequent calculation of 
#'  the statistics. The exception to this design are the correlation structure 
#'  based statistics as they require all entries of the correlation structure 
#'  matrices for that module. As modules can be quite large (e.g., 10,000 
#'  variables), this would be dramatically increase memory overhead to apply 
#'  this strategy.
#' 
#' @references
#'  \enumerate{
#'     \item{
#'       Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my
#'       network module preserved and reproducible?} PLoS Comput. Biol. 
#'       \strong{7}, e1001057 (2011). 
#'     }
#'  }
#'  
#' @return 
#' a vector containing the values for all the module preservation
#' statistics: \emph{avg.weight}, \emph{cor.degree}, 
#' \emph{coherence}, \emph{avg.contrib}, \emph{cor.contrib}, \emph{cor.cor}, and 
#' \emph{avg.cor}.
calcStats <- function(
  discProps, testProps, discCor, discIndices, testCor, testIndices
) {
  cList <- corStats(discCor, discIndices, testCor, testIndices)
  stats <- c(
    avg.weight = testProps[["averageEdgeWeight"]],
    cor.degree = cor(discProps[["weightedDegree"]], testProps[["weightedDegree"]]),
    cor.cor = cor(cList[["cor.discovery"]], cList[["cor.test"]]),
    avg.cor = cList[["corDensity"]]
  )
  if ("moduleCoherence" %in% names(testProps)) { # Detect if data has been provided
    stats <- c(
      stats, coherence = testProps[["moduleCoherence"]],
      avg.contrib = mean(sign(discProps[["nodeContribution"]]) * testProps[["nodeContribution"]]),
      cor.contrib = cor(discProps[["nodeContribution"]], testProps[["nodeContribution"]])
    )
  }
  stats
}

#' Calculate a contigency table of module overlap between datasets
#' 
#' @param moduleAssignments processed user input for the 'moduleAssignments'
#'  argument (see \link{processInput})
#' @param modules user processed input for the 'modules' argument.
#' @param network user processd input for the 'network' argument.
#' @param di the discovery dataset to use
#' @param ti the test dataset to use
#' 
#' @return a list containing a contigency table, 
contingencyTable <- function(moduleAssignments, modules, network, di, ti) {
  # To simplify later function calls, we need to get a vector of module
  # assignments only for (a) modules of interest and (b) the variables
  # present in both datasets for those modules.
  overlapVars <- intersect(
    colnames(network[[di]]), 
    colnames(network[[ti]])
  )
  overlapAssignments <- moduleAssignments[[di]][overlapVars]
  
  overlapAssignments <- overlapAssignments %sub_in% modules[[di]] 
  overlapModules <- unique(overlapAssignments)
  overlapModules <- overlapModules[orderAsNumeric(overlapModules)]
  
  if (length(overlapAssignments) == 0) {
    warning(
      "No variables composing the modules of interest are present in", 
      " the test dataset"
    )
    next
  }
  
  # How many variables are present in the test dataset for the modules 
  # of interest?
  varsPres <- table(overlapAssignments)
  modulesWithNoOverlap <- modules[[di]] %sub_nin% overlapAssignments
  varsPres <- c(varsPres, rep(0, length(modulesWithNoOverlap)))
  names(varsPres)[names(varsPres) == ""] <- modulesWithNoOverlap
  varsPres <- varsPres[orderAsNumeric(names(varsPres))]
  
  # What proportion?
  moduleSizes <- table(moduleAssignments[[di]])
  moduleSizes <- moduleSizes[names(moduleSizes) %sub_in% modules[[di]]]
  propVarsPres <- varsPres / moduleSizes
  propVarsPres <- propVarsPres[orderAsNumeric(names(propVarsPres))]
  
  # Calculate some basic cross-tabulation statistics so we can assess 
  # which modules in both datasets map to each other, if module
  # detection has also been performed for the test network
  contingency <- NULL
  if (!is.null(moduleAssignments[[ti]])) {
    # Get total number of nodes from each discovery subset in each test subset 
    contingency <- table(
      moduleAssignments[[di]][overlapVars], 
      moduleAssignments[[ti]][overlapVars]
    )
    # filter on subsets the user cares about
    contingency <- contingency[modules[[di]],,drop=FALSE]
    
    # Order numerically if relevant
    contingency <- contingency[
      orderAsNumeric(rownames(contingency)),
      orderAsNumeric(colnames(contingency)),
      drop=FALSE
    ]
    
    # add in the module sizes from the respective datasets
    contingency <- cbind(rowSums(contingency), contingency)
    testSizes <- table(moduleAssignments[[ti]][overlapVars])
    testSizes <- testSizes[colnames(contingency)]
    
    contingency <- rbind(
      testSizes[colnames(contingency)], 
      contingency
    )
    rownames(contingency)[1] <- "size"
    colnames(contingency)[1] <- "size"
  }
  
  return(list(
    contingency=contingency, propVarsPres=propVarsPres, overlapVars=overlapVars,
    varsPres=varsPres, overlapModules=overlapModules,
    overlapAssignments=overlapAssignments
  ))
}
