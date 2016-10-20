#' Calculate a contigency table of module overlap between datasets
#' 
#' @param modAssignments a list where the first element is the 
#'  'moduleAssignments' vector in the discovery dataset, and the second element
#'  is the 'moduleAssignments vector in the test dataset.
#' @param mods the 'modules' vector for the discovery dataset.
#' @param tiNodelist a vector of node IDs in the test dataset.
#' 
#' @return a list containing a contigency table, 
contingencyTable <- function(modAssignments, mods, tiNodelist) {
  # To simplify later function calls, we need to get a vector of module
  # assignments only for (a) modules of interest and (b) the variables
  # present in both datasets for those modules.
  overlapVars <- intersect(names(modAssignments[[1]]), tiNodelist)
  overlapAssignments <- modAssignments[[1]][overlapVars]
  
  overlapAssignments <- overlapAssignments %sub_in% mods
  overlapModules <- unique(overlapAssignments)
  overlapModules <- overlapModules[orderAsNumeric(overlapModules)]
  
  # How many variables are present in the test dataset for the modules 
  # of interest?
  varsPres <- table(overlapAssignments)
  
  modulesWithNoOverlap <- mods %sub_nin% overlapAssignments
  varsPres <- c(varsPres, rep(0, length(modulesWithNoOverlap)))
  names(varsPres)[names(varsPres) == ""] <- modulesWithNoOverlap
  varsPres <- varsPres[orderAsNumeric(names(varsPres))]
  
  if (any(varsPres == 0)) {
    noNodes <- names(varsPres[varsPres == 0])
    warning(
      "None of the nodes in module(s) ", 
      paste(paste0('"', noNodes, '"'), collapse=", "),
      " are present in the test dataset",
      immediate. = TRUE
    )
  }
  
  # What proportion?
  moduleSizes <- table(modAssignments[[1]])
  moduleSizes <- moduleSizes[names(moduleSizes) %sub_in% mods]
  propVarsPres <- varsPres / moduleSizes[names(varsPres)]

  # Calculate some basic cross-tabulation statistics so we can assess 
  # which modules in both datasets map to each other, if module
  # detection has also been performed for the test network
  contingency <- NULL
  if (!is.null(modAssignments[[2]])) {
    # Get total number of nodes from each discovery subset in each test subset 
    contingency <- table(
      modAssignments[[1]][overlapVars], 
      modAssignments[[2]][overlapVars]
    )
    # filter on subsets the user cares about
    contingency <- contingency[mods,,drop=FALSE]
    
    # Order numerically if relevant
    contingency <- contingency[
      orderAsNumeric(rownames(contingency)),
      orderAsNumeric(colnames(contingency)),
      drop=FALSE
    ]
    
    # add in the module sizes from the respective datasets
    contingency <- cbind(rowSums(contingency), contingency)
    testSizes <- table(modAssignments[[2]][overlapVars])
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
