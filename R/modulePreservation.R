#' Replication and preservation of network modules across datasets
#' 
#' Quantify the preservation of network modules (sub-graphs) in an independent
#' dataset through permutation testing on module topology. Seven network
#' statistics (see details) are calculated for each module and then tested by
#' comparing to distributions generated from their calculation on random subsets
#' in the test dataset.
#'
#' @inheritParams common_params
#' @inheritParams simplify_param
#' 
#' @param selfPreservation logical; if \code{FALSE} (default) then module 
#'  preservation analysis will not be performed within a dataset (\emph{i.e.} 
#'  where the \code{discovery} and \code{test} datasets are the same).
#' @param nThreads number of threads to parallelise the calculation of network 
#'   properties over. Automatically determined as the number of cores - 1 if 
#'   not specified.
#' @param nPerm number of permutations to use. If not specified, the number of 
#'  permutations will be automatically determined (see details).
#' @param null variables to include when generating the null distributions. 
#'  Must be either "overlap" or "all" (see details).
#' @param alternative The type of module preservation test to perform. Must be 
#'   one of "greater" (default), "less" or "two.sided" (see details).
#'  
#' @details
#'  \subsection{Input data structures:}{
#'   The preservation of network modules in a second dataset is quantified by
#'   measuring the preservation of topological properties between the
#'   \emph{discovery} and \emph{test} datasets. These properties are calculated
#'   not only from the interaction networks inferred in each dataset, but also
#'   from the data used to infer those networks (e.g. gene expression data) as
#'   well as the correlation structure between variables/nodes. Thus, all
#'   functions in the \code{NetRep} package have the following arguments:
#'   \itemize{
#'     \item{\code{network}:}{
#'       a list of interaction networks, one for each dataset.
#'     }
#'     \item{\code{data}:}{
#'       a list of data matrices used to infer those networks, one for each 
#'       dataset.
#'     }
#'     \item{\code{correlation}:}{
#'      a list of matrices containing the pairwise correlation coefficients 
#'      between variables/nodes in each dataset.
#'     } 
#'     \item{\code{moduleAssignments}:}{
#'      a list of vectors, one for each \emph{discovery} dataset, containing 
#'      the module assignments for each node in that dataset.
#'     }
#'     \item{\code{modules}:}{
#'      a list of vectors, one for each \emph{discovery} dataset, containing
#'      the names of the modules from that dataset to analyse.  
#'     }
#'     \item{\code{discovery}:}{
#'       a vector indicating the names or indices of the previous arguments' 
#'       lists to use as the \emph{discovery} dataset(s) for the analyses.
#'     }
#'     \item{\code{test}:}{
#'       a list of vectors, one vector for each \emph{discovery} dataset, 
#'       containing the names or indices of the \code{network}, \code{data}, and 
#'       \code{correlation} argument lists to use as the \emph{test} dataset(s) 
#'       for the analysis of each \emph{discovery} dataset.
#'     }
#'   }
#'   
#'   The formatting of these arguments is not strict: each function will attempt
#'   to make sense of the user input. For example, if there is only one 
#'   \code{discovery} dataset, then input to the \code{moduleAssigments} and 
#'   \code{test} arguments may be vectors, rather than lists. 
#' }
#' \subsection{Analysing large datasets:}{
#'   Matrices in the \code{network}, \code{data}, and \code{correlation} lists
#'   can be supplied as \code{\link{disk.matrix}} objects. This class allows 
#'   matrix data to be kept on disk and loaded as required by \pkg{NetRep}. 
#'   This dramatically decreases memory usage: the matrices for only one 
#'   dataset will be kept in RAM at any point in time.
#'   
#'   Additional memory usage of the permutation procedure is directly
#'   proportional to the sum of module sizes squared multiplied by the number 
#'   of threads. Very large modules may result in significant additional memory
#'   usage per core due to extraction of the correlation coefficient sub-matrix
#'   at each permutation.
#' }
#' \subsection{Module Preservation Statistics:}{
#'  Module preservation is assessed through seven module preservation statistics,
#'  each of which captures a different aspect of a module's topology; \emph{i.e.}
#'  the structure of the relationships between its nodes \emph{(1,2)}. Below is
#'  a description of each statistic, what they seek to measure, and where their
#'  interpretation may be inappropriate. 
#'  
#'  The \emph{module coherence} (\code{'coherence'}), \emph{average node 
#'  contribution} (\code{'avg.contrib'}), and \emph{concordance of node 
#'  contribution} (\code{'cor.contrib'}) are all calculated from the data used 
#'  to infer the network (provided in the \code{'data'} argument). They are 
#'  calculated from the module's \emph{summary profile}. This is the eigenvector
#'  of the 1st principal component across all observations for every node
#'  composing the module. For gene coexpression modules this can be interpreted
#'  as a "summary expression profile". It is typically referred to as the
#'  "module eigengene" in the weighted gene coexpression network analysis
#'  literature \emph{(4)}.
#'  
#'  The \emph{module coherence} (\code{'coherence'}) quantifies the proportion 
#'  of module variance explained by the module's "summary profile". The higher
#'  this value, the more "coherent" the data is, \emph{i.e.} the more similar
#'  the observations are nodes for each sample. With the default alternate
#'  hypothesis, a small permutation \emph{P}-value indicates that the module is
#'  more coherent than expected by chance.
#'  
#'  The \emph{average node contribution} (\code{'avg.contrib'}) and 
#'  \emph{concordance of node contribution} (\code{'cor.contrib'}) are calculated 
#'  from the \emph{node contribution}, which quantifies how similar each node is 
#'  to the modules's \emph{summary profile}. It is calculated as the Pearson
#'  correlation coefficient between each node and the module summary profile. In
#'  the weighted gene coexpression network literature it is typically called the
#'  "module membership" \emph{(2)}.
#'  
#'  The \emph{average node contribution} (\code{'avg.contrib'}) quantifies how
#'  similar nodes are to the module summary profile in the test dataset. Nodes
#'  detract from this score where the sign of their node contribution flips 
#'  between the discovery and test datasets, \emph{e.g.} in the case of 
#'  differential gene expression across conditions. A high \emph{average node
#'  contribution} with a small permutation \emph{P}-value indicates that the
#'  module remains coherent in the test dataset, and that the nodes are acting
#'  together in a similar way.  
#'  
#'  The \emph{concordance of node contribution} (\code{'cor.contrib'}) measures 
#'  whether the relative rank of nodes (in terms of their node contribution) is 
#'  preserved across datasets. If a module is coherent enough that all nodes 
#'  contribute strongly, then this statistic will not be meaningful as its value
#'  will be heavily influenced by tiny variations in node rank. This can be
#'  assessed through visualisation of the module topology (see 
#'  \code{\link{plotContribution}}.) Similarly, a strong
#'  \code{'cor.contrib'} is unlikely to be meaningful if the
#'  \code{'avg.contrib'} is not significant.
#'  
#'  The \emph{concordance of correlation strucutre} (\code{'cor.cor'}) and 
#'  \emph{density of correlation structure} (\code{'avg.cor'}) are calculated 
#'  from the user-provided correlation structure between nodes (provided in the 
#'  \code{'correlation'} argument). This is referred to as "coexpression" when
#'  calculated on gene expression data.
#'  
#'  The \code{'avg.cor'} measures how strongly nodes within a module are 
#'  correlation on average in the test dataset. This average depends on the 
#'  correlation coefficients in the discovery dataset: the score is penalised 
#'  where correlation coefficients change in sign between datasets. A high 
#'  \code{'avg.cor'} with a small permutation \emph{P}-value indicates that the 
#'  module is (a) more strongly correlated than expected by chance for a module 
#'  of the same size, and (b) more consistently correlated with respect to the 
#'  discovery dataset than expected by chance.
#'  
#'  The \code{'cor.cor'} measures how similar the correlation coefficients are 
#'  across the two datasets. A high \code{'cor.cor'} with a small permutation 
#'  \emph{P}-value indicates that the correlation structure within a module is 
#'  more similar across datasets than expected by chance. If all nodes within a 
#'  module are very similarly correlated then this statistic will not be 
#'  meaningful, as its value will be heavily influenced by tiny, non-meaningful, 
#'  variations in correlation strength. This can be assessed through
#'  visualisation of the module topology (see \code{\link{plotCorrelation}}.)
#'  Similarly, a strong \code{'cor.cor'} is unlikely to be meaningful if the
#'  \code{'avg.cor'} is not significant.
#'  
#'  The \emph{average edge weight} (\code{'avg.weight'}) and \emph{concordance
#'  of weighted degree} (\code{'cor.degree'}) are both calculated from the 
#'  interaction network (provided as adjacency matrices to the \code{'network'}
#'  argument). 
#'  
#'  The \code{'avg.weight'} measures the average connection strength between 
#'  nodes in the test dataset. In the weighted gene coexpression network 
#'  literature this is typically called the "module density" \emph{(2)}. A high
#'  \code{'avg.weight'} with a small permutation \emph{P}-value indicates that
#'  the module is more strongly connected in the test dataset than expected by
#'  chance. 
#'  
#'  The \code{'cor.degree'} calculates whether the relative rank of each node's 
#'  \emph{weighted degree} is similar across datasets. The \emph{weighted
#'  degree} is calculated as the sum of a node's edge weights to all other nodes
#'  in the module. In the weighted gene coexpression network literature this is 
#'  typically called the "intramodular connectivity" \emph{(2)}. This statistic 
#'  will not be meaningful where all nodes are connected to each other with 
#'  similar strength, as its value will be heavily influenced by tiny,
#'  non-meaningful, variations in weighted degree. This can be assessed through
#'  visualisation of the module topology (see \code{\link{plotDegree}}.)
#'  
#'  Both the \code{'avg.weight'} and \code{'cor.degree'} assume edges are 
#'  weighted, and that the network is densely connected. Note that for sparse 
#'  networks, edges with zero weight are included when calculating both
#'  statistics. Only the magnitude of the weights, not their sign, contribute to
#'  the score. If the network is \emph{unweighted}, \emph{i.e.} edges indicate
#'  presence or absence of a relationship, then the \code{'avg.weight'} will be
#'  the proportion of the number of edges to the total number of possible edges
#'  while the \emph{weighted degree} simply becomes the \emph{degree}. A high
#'  \code{'avg.weight'} in this case measures how interconnected a module is in
#'  the test dataset. A high \emph{degree} indicates that a node is connected to
#'  many other nodes. The interpretation of the \code{'cor.degree'} remains
#'  unchanged between weighted and unweighted networks. If the network is
#'  directed the interpretation of the \code{'avg.weight'} remains unchanged,
#'  while the \emph{cor.degree} will measure the concordance of the node
#'  \emph{in-}degree in the test network. To measure the \emph{out-}degree
#'  instead, the adjacency matrices provided to the \code{'network'} argument
#'  should be transposed.
#' }
#' \subsection{Sparse data:}{
#'  Caution should be used when running \code{NetRep}
#'  on sparse data (\emph{i.e.} where there are many zero values in the data 
#'  used to infer the network). For this data, the \emph{average node contribution} 
#'  (\code{'avg.contrib'}), \emph{concordance of node contribution} 
#'  (\code{'cor.contrib'}), and \emph{module coherence} (\code{'coherence'})
#'  will all be systematically underestimated due to their reliance on the 
#'  Pearson correlation coefficient to calculate the \emph{node contribution}.
#'  
#'  Care should also be taken to use appropriate methods for inferring the
#'  correlation structure when the data is sparse for the same reason.
#' }
#' \subsection{Proportional data:}{
#'  Caution should be used when running \code{NetRep} on proportional data (
#'  \emph{i.e.} where observations across samples all sum to the same value, 
#'  \emph{e.g.} 1). For this data, the \emph{average node contribution} 
#'  (\code{'avg.contrib'}), \emph{concordance of node contribution} 
#'  (\code{'cor.contrib'}), and \emph{module coherence} (\code{'coherence'})
#'  will all be systematically overestimated due to their reliance on the 
#'  Pearson correlation coefficient to calculate the \emph{node contribution}.
#'  
#'  Care should also be taken to use appropriate methods for inferring the
#'  correlation structure from proportional data for the same reason.
#' }
#' \subsection{Hypothesis testing:}{
#'  Three alternative hypotheses are available. "greater", the default, tests
#'  whether each module preservation statistic is larger than expected by 
#'  chance. "lesser" tests whether each module preservation statistic is smaller
#'  than expected by chance, which may be useful for identifying modules that
#'  are extremely different in the \emph{test} dataset. "two.sided" can be used
#'  to test both alternate hypotheses.
#'  
#'  To determine whether a module preservation statistic deviates from chance, a
#'  permutation procedure is employed. Each statistic is calculated between the
#'  module in the \emph{discovery} dataset and \code{nPerm} random subsets of
#'  the same size in the \emph{test} dataset in order to assess the distribution
#'  of each statistic under the null hypothesis. 
#'  
#'  Two models for the null hypothesis are available: "overlap", the default, 
#'  only nodes that are present in both the \emph{discovery} and \emph{test}
#'  networks are used when generating null distributions. This is appropriate
#'  under an assumption that nodes that are present in the \emph{test} dataset, 
#'  but not present in the \emph{discovery} dataset, are unobserved: that is,
#'  they may fall in the module(s) of interest in the \emph{discovery} dataset
#'  if they were to be measured there. Alternatively, "all" will use all nodes
#'  in the \emph{test} network when generating the null distributions.
#'   
#'  The number of permutations required for any given significance threshold is 
#'  approximately 1 / the desired significance for one sided tests, and double 
#'  that for two-sided tests. This can be calculated with 
#'  \code{\link{requiredPerms}}. When \code{nPerm} is not specified, the number 
#'  of permutations is automatically calculated as the number required for a 
#'  Bonferroni corrected significance threshold adjusting for the total number 
#'  of tests for each statistic, i.e. the total number of modules to be analysed
#'  multiplied by the number of \emph{test} datasets each module is tested in. 
#'  Although assessing the replication of a small numberof modules calls for 
#'  very few permutations, we recommend using no fewer than 1,000 as fewer 
#'  permutations are unlikely to generate representative null distributions. 
#'  \strong{Note:} the assumption used by \code{\link{requiredPerms}} to 
#'  determine the correct number of permtutations breaks down when assessing the
#'  preservation of modules in a very small dataset (e.g. gene sets in a dataset
#'  with less than 100 genes total). However, the reported p-values will still
#'  be accurate (see \code{\link{permutationTest}}) \emph{(3)}.
#' }
#' 
#' @references 
#'   \enumerate{
#'     \item{
#'      Ritchie, S.C., \emph{et al.}, \emph{A scalable permutation approach
#'      reveals replication and preservation patterns of gene coexpression
#'      modules}. bioRxiv. 029553 (2015).
#'     }
#'     \item{
#'       Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my
#'       network module preserved and reproducible?} PLoS Comput. Biol. 
#'       \strong{7}, e1001057 (2011). 
#'     }
#'     \item{
#'       Phipson, B. & Smyth, G. K. \emph{Permutation P-values should never be 
#'       zero: calculating exact P-values when permutations are randomly drawn.}
#'       Stat. Appl. Genet. Mol. Biol. \strong{9}, Article39 (2010). 
#'     }
#'     \item{
#'       Langfelder, P. & Horvath, S. \emph{WGCNA: an R package for weighted 
#'       correlation network analysis.} BMC Bioinformatics \strong{9}, 559 
#'       (2008).
#'     }
#'   }
#'   
#' @return
#'  A nested list structure. At the top level, the list has one element per 
#'  \code{'discovery'} dataset. Each of these elements is a list that has one
#'  element per \code{'test'} dataset analysed for that \code{'discovery'} 
#'  dataset. Each of these elements is also a list, containing the following
#'  objects:
#'  \itemize{
#'    \item{\code{observed}:}{
#'      A matrix of the observed values for the module preservation statistics.
#'      Rows correspond to modules, and columns to the module preservation
#'      statistics.
#'    }
#'    \item{\code{nulls}:}{
#'      A three dimensional array containing the values of the module 
#'      preservation statistics evaluated on random permutation of module 
#'      assignment in the test network. Rows correspond to modules, columns to
#'      the module preservation statistics, and the third dimension to the 
#'      permutations.
#'    }
#'    \item{\code{p.values}:}{
#'      A matrix of p-values for the \code{observed} module preservation 
#'      statistics as evaluated through a permutation test using the 
#'      corresponding values in \code{nulls}.
#'    }
#'    \item{\code{nVarsPresent}:}{
#'      A vector containing the number of variables that are present in the test
#'      dataset for each module.
#'    }
#'    \item{\code{propVarsPresent}:}{
#'      A vector containing the proportion of variables present in the test dataset
#'      for each module. Modules where this is less than 1 should be 
#'      investigated further before making judgements about preservation to 
#'      ensure that the missing variables are not the most connected ones.
#'    }
#'    \item{\code{contingency}:}{ 
#'      If \code{moduleAssignments} are present for both the \emph{discovery}
#'      and \emph{test} datasets, then a contingency table showing the overlap
#'      between modules across datasets is returned. Rows correspond to modules
#'      in the \emph{discovery} dataset, columns to modules in the \emph{test}
#'      dataset.
#'    }
#'  }
#'  When \code{simplify = TRUE} then the simplest possible structure will be 
#'  returned. E.g. if module preservation is tested in only one dataset, then
#'  the returned list will have only the above elements. 
#'  
#'  When \code{simplify = FALSE} then a nested list of datasets will always be 
#'  returned, i.e. each element at the top level and second level correspond to
#'  a dataset, e.g. \code{results[["Dataset1"]][["Dataset2"]]} indicates an 
#'  analysis where modules discovered in "Dataset1" are assessed for 
#'  preservation in "Dataset2". Dataset comparisons which have not been 
#'  assessed will contain \code{NULL}.
#'  
#' @seealso 
#'   Functions for: 
#'   \link[=plotModule]{visualising network modules},
#'   \link[=networkProperties]{calculating module topology}, 
#'   \link[=permutationTest]{calculating permutation test P-values}, and 
#'   \link[=combineAnalyses]{splitting computation over multiple machines}.
#'   
#' @examples
#' # load in example data, correlation, and network matrices for a discovery and test dataset:
#' data("NetRep")
#' 
#' # Set up input lists for each input matrix type across datasets. The list
#' # elements can have any names, so long as they are consistent between the
#' # inputs.
#' network_list <- list(discovery=discovery_network, test=test_network)
#' data_list <- list(discovery=discovery_data, test=test_data)
#' correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
#' labels_list <- list(discovery=module_labels)
#' 
#' # Assess module preservation.
#' preservation <- modulePreservation(
#'  network=network_list, data=data_list, correlation=correlation_list, 
#'  moduleAssignments=labels_list, nPerm=10000, discovery="discovery", 
#'  test="test", nThreads=2
#' )
#' 
#' @importFrom parallel detectCores
#' @import foreach
#' @import RhpcBLASctl
#' @export
modulePreservation <- function(
  network, data, correlation, moduleAssignments, modules=NULL, 
  backgroundLabel="0", discovery=1, test=2, selfPreservation=FALSE,
  nThreads=NULL, nPerm=NULL, null="overlap", alternative="greater", 
  simplify=TRUE, verbose=TRUE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  # Identify the null hypothesis to use 
  nullModels <- c("overlap", "all")
  if (is.na(pmatch(null, nullModels))) {
    stop("overlap must match one of ", paste(nullModels, collapse=" "))
  }
  model <- nullModels[pmatch(null, nullModels)]
  
  # Identify the alternate hypothesis to use
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  alternative <- validAlts[altMatch]
  
  # Validate 'nPerm'. If 'NULL', we need to process the rest of the input to
  # determine.
  if (!is.null(nPerm) && (!is.numeric(nPerm) || length(nPerm) > 1 || nPerm < 1)) {
    stop("'nPerm' must be a single number > 1")
  }
  
  # Validate 'nThreads'
  maxThreads <- detectCores()
  if (is.null(nThreads)) {
    if (is.na(maxThreads)) {
      stop("'nThreads' must always be supplied by the user on this machine")
    }
    nThreads <- maxThreads - 1; # Leave a core for interactive use
  }
  
  if (!is.numeric(nThreads) || length(nThreads) > 1 || nThreads < 1)
    stop("'nThreads' must be a single number greater than 0")
  
  if (!is.na(maxThreads) && nThreads > maxThreads) {
    stop(    
      "Number of threads requested (", nThreads, ") exceeds the reported ",
      "maximum number of concurrent threads supported by the current ", 
      "hardware (", maxThreads, ")."
    )
  }
  
  # Disable implicit parallelism (i.e. through multithreaded BLAS libraries).
  # The permutation procedure will spawn its own threads.
  oldOMPThreads <- omp_get_max_threads()
  oldBLASThreads <- blas_get_num_procs()
  
  omp_set_num_threads(1)
  blas_set_num_threads(1)
  
  # Restore to previous state
  on.exit({
    omp_set_num_threads(oldOMPThreads)
    blas_set_num_threads(oldBLASThreads)
  }, add=TRUE)

  # For the other functions this isn't necessary.
  if (missing(moduleAssignments))
    moduleAssignments <- NULL
  if (is.null(moduleAssignments))
    stop("'moduleAssignments' must be provided for each 'discovery' dataset")
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose)
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  discovery <- finput$discovery
  test <- finput$test
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  nodelist <- finput$nodelist

  # If NULL, automatically determine.
  if (is.null(nPerm)) {
    # If missing, set as the required number for Bonferroni correction.
    # Bonferonni correct for the total number of modules, multiplied the number
    # of datasets each module is tested in.
    multiplier <- sum(sapply(modules, length) * sapply(test, length))
    nPerm <- max(10000, requiredPerms(0.05/multiplier))
  }

  vCat(verbose, 0, "Input ok!")
  
  # Set up return list
  res <- foreach(di = seq_len(nDatasets)) %do% {
    res2 <- foreach(ti = seq_len(nDatasets)) %do% {}
    names(res2) <- datasetNames
    return(res2)
  } 
  names(res) <- datasetNames
  
  #-----------------------------------------------------------------------------
  # Set up variables for running module preservation analysis
  #-----------------------------------------------------------------------------
  
  # The following declarations are for iterators declared inside each foreach 
  # loop. Declarations are required to satisfy NOTES generated by R CMD check, 
  # and also serve as useful documentation for the reader of the source code.
  di <- NULL # discovery dataset 
  ti <- NULL # test dataset 
  mi <- NULL # iterator over the modules
  si <- NULL # iterator over the statistics
  pi <- NULL # iterator over the permutations
  
  #-----------------------------------------------------------------------------
  # Pairwise iteration over the datasets
  #-----------------------------------------------------------------------------
  for (di in discovery) {
    for (ti in test[[di]]) {
      if (!selfPreservation && di == ti) {
        vCat(
          verbose, 0, sep="", "skipping module preservation analysis for modules",
          " from dataset ", '"', datasetNames[di], '"', " within dataset ", '"', 
          datasetNames[di], '"', "."
        )
        next
      }
      tryCatch({
        vCat(
          verbose, 0, sep="", 
          "Calculating preservation of network subsets from dataset ", '"', 
          datasetNames[di], '"', " in dataset ", '"', datasetNames[ti], '"', "."
        )
        #----------------------------------------------------------------------
        # Calculate the overlap between datasets
        #----------------------------------------------------------------------
        ct <- contingencyTable(moduleAssignments[c(di, ti)], modules[[di]], 
                               nodelist[[ti]])
        contingency <- ct$contingency
        propVarsPres <- ct$propVarsPres
        overlapVars <- ct$overlapVars
        varsPres <- ct$varsPres
        overlapModules <- ct$overlapModules
        overlapAssignments <- ct$overlapAssignments
        
        nStatistics <- ifelse(!is.null(data[[di]]) && !is.null(data[[ti]]), 7, 4)
        nModules <- length(overlapModules)
        
        #----------------------------------------------------------------------
        # Calculate the intermediate properties of the discovery dataset
        #----------------------------------------------------------------------
        # These are needed at every permutation, but we can cut runtime by 
        # calculating them once, and cut memory by loading and unloading the
        # the discovery dataset if provided as 'disk.matrix' objects.
        anyDM <- any.disk.matrix(data[[di]], correlation[[di]], network[[di]])
        vCat(verbose && anyDM, 1, 'Loading matrices of dataset "', 
             datasetNames[di], '" into RAM...', sep="")
        if (!is.null(data[[di]]) && !is.null(data[[ti]])) {
          discovery_data <- loadIntoRAM(data[[di]])
        } else {
          discovery_data <- NULL
        }
        discovery_correlation <- loadIntoRAM(correlation[[di]])
        discovery_network <- loadIntoRAM(network[[di]])

        # Calculate the intermediate properties
        vCat(verbose, 1, 'Pre-computing intermediate properties in dataset "',
             datasetNames[di], '"...', sep="")
        if (is.null(data[[di]]) || is.null(data[[ti]])) {
          discProps <- IntermediatePropertiesNoData(
            discovery_correlation, discovery_network, nodelist[[ti]],
            moduleAssignments[[di]], modules[[di]]
          )
        } else {
          discProps <- IntermediateProperties(
            discovery_data, discovery_correlation, discovery_network,
            nodelist[[ti]], moduleAssignments[[di]], modules[[di]]
          )
        }

        # Free up memory
        vCat(verbose && anyDM, 1, "Unloading matrices...")
        rm(discovery_data, discovery_correlation, discovery_network)
        gc()
        
        #----------------------------------------------------------------------
        # Run the permutation procedure
        #----------------------------------------------------------------------
        # Load matrices into RAM if they are 'disk.matrix' objects.
        anyDM <- any.disk.matrix(data[[ti]], correlation[[ti]], network[[ti]])
        vCat(verbose && anyDM, 1, 'Loading matrices of dataset "', 
             datasetNames[ti], '" into RAM...', sep="")
        if (!is.null(data[[di]]) && !is.null(data[[ti]])) {
          test_data <- loadIntoRAM(data[[ti]])
        } else {
          test_data <- NULL
        }
        test_correlation <- loadIntoRAM(correlation[[ti]])
        test_network <- loadIntoRAM(network[[ti]])

        # Run the permutation procedure
        if (is.null(data[[di]]) || is.null(data[[ti]])) {
          perms <- PermutationProcedureNoData(
            discProps, test_correlation, test_network, moduleAssignments[[di]], 
            modules[[di]], nPerm, nThreads, model, verbose, vCat
          )
        } else {
          perms <- PermutationProcedure(
            discProps, test_data, test_correlation, test_network, 
            moduleAssignments[[di]], modules[[di]], nPerm, nThreads, model, 
            verbose, vCat
          )
        }
        observed <- perms$observed
        nulls <- perms$nulls
        
        vCat(verbose && anyDM, 1, "Unloading matrices...")
        # Free up memory
        rm(test_data, test_correlation, test_network)
        gc()
        
        #---------------------------------------------------------------------
        # Calculate permutation p-value
        #---------------------------------------------------------------------
        vCat(verbose, 1, "Calculating P-values...")
        if (model == 'overlap') {
          totalSize <- length(overlapVars)
        } else {
          totalSize <- ncol(correlation[[ti]])
        }
        p.values <- permutationTest(nulls, observed, varsPres, totalSize, alternative)
        
        #---------------------------------------------------------------------
        # Collate results
        #---------------------------------------------------------------------
        vCat(verbose, 1, "Collating results...")
    
        res[[di]][[ti]] <- list(
          nulls = nulls,
          observed = observed,
          p.values = p.values,
          nVarsPresent = varsPres,
          propVarsPresent = propVarsPres,
          totalSize = totalSize,
          alternative = alternative,
          contingency = contingency
        )
        # remove NULL outputs
        res[[di]][[ti]] <- res[[di]][[ti]][
          !sapply(res[[di]][[ti]], is.null)
        ]
        
        gc()
      }, error=function(e) {
        warning(
          "Failed with error:\n", e$message, "\nSkipping to next comparison",
          immediate. = TRUE
        )
      })
    }
  }
  
  # Simplify the output data structure where possible
  if (simplify) {
    res <- simplifyList(res, depth=2)
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
  res
}
