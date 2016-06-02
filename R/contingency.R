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
  
  # How many variables are present in the test dataset for the modules 
  # of interest?
  varsPres <- table(overlapAssignments)
  
  modulesWithNoOverlap <- modules[[di]] %sub_nin% overlapAssignments
  varsPres <- c(varsPres, rep(0, length(modulesWithNoOverlap)))
  names(varsPres)[names(varsPres) == ""] <- modulesWithNoOverlap
  varsPres <- varsPres[orderAsNumeric(names(varsPres))]
  
  if (any(varsPres == 0)) {
    noNodes <- names(varsPres[varsPres == 0])
    warning(
      "None of the nodes in module(s) ", 
      paste(paste0('"', noNodes, '"'), collapse=", "),
      " are present in the test dataset"
    )
  }
  
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
