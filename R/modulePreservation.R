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
#'  preservation analysis will not be performed where the \code{discovery} and
#'  \code{test} datasets are the same.
#' @param nPerm number of permutations to use. If not specified, the number of 
#'  permutations will be automatically determined (see details).
#' @param nCores number of cores to parallelise the permutation procedure over.
#'  Ignored if the user has already registered a parallel backend. If 
#'  \code{NULL} (default) the maximum number of cores on the machine will be 
#'  used.
#' @param null variables to include when generating the null distributions. 
#'  Must be either "overlap" or "all" (see details).
#' @param alternative The type of module preservation test to perform. Must be 
#'   one of "greater" (default), "less" or "two.sided" (see details).
#' @param statCorMethod character vector indicating method to use when calculating 
#'   the correlation based statistics (see details). Must be one of "pearson", 
#'   "spearman", or "kendall". If the WGCNA package is installed then "bicor" 
#'   may also be specified as an option (see \code{\link[WGCNA]{bicor}}). 
#'  
#' @details
#'  \subsection{Input data structure:}{
#'   This function allows for input data formatted in a number of ways. Where 
#'   there are multiple datasets of interest (e.g. multiple tissues, locations, 
#'   or a discovery dataset and an independent test dataset) the arguments 
#'   \code{data}, \code{correlation}, and \code{network} should be
#'   \code{\link[=list]{lists}} where each element contains the matrix data for 
#'   each respective dataset. Alternatively, if only one dataset is of interest, 
#'   the \code{data}, \code{correlation}, and \code{network} arguments
#'   will also each accept a 'matrix' object.
#'   
#'   Similarly, the \code{moduleAssignments} argument expects a list of named
#'   vectors, which denote the module each variable belongs to in the discovery
#'   dataset. If module discovery has only been performed in one dataset, then 
#'   the \code{moduleAssignments} argument will also accept a named vector.
#'   
#'   The \code{discovery} arguments specifies which dataset the \code{modules} 
#'   of interest were discovered in, and the \code{test} argument specifies 
#'   which dataset to calculate the network properties in. These arguments are
#'   ignored if data is provided for only one dataset.
#' }
#' \subsection{'bigMatrix' vs. 'matrix' input data:}{
#'   Although the function expects \code{\link[=bigMatrix-class]{bigMatrix}} 
#'   data, regular 'matrix' objects are also accepted. In this case, the 
#'   'matrix' data is temporarily converted to 'bigMatrix' by the function. This
#'   conversion process involves writing out each matrix as a binary file on 
#'   disk, which can take a long time for large datasets. It is strongly 
#'   recommended for the user to store their data as 'bigMatrix' objects, as the
#'   \link{modulePreservation} function, \link[=plotModule]{plotting} 
#'   \link[=plotTopology]{functions}, \link[=nodeOrder]{node} and 
#'   \link[=sampleOrder]{sample} ordering functions also expect 'bigMatrix'
#'   objects. Further, 'bigMatrix' objects have a number of benefits, including 
#'   instantaneous load time from any future R session, and parallel access from
#'   mutliple independent R sessions. Methods are provided for 
#'   converting to, loading in, and writing out \code{\link{bigMatrix}} objects.
#' }
#' \subsection{Memory usage:}{
#'   A trade off has been made between memory usage and computation time. 
#'   'modulePreservation' has a large memory overhead as it requires 
#'   pre-computed correlation and network matrices for each dataset. However,
#'   these are stored in shared memory, which means that each parallel process
#'   can independently access this memory. There is very little memory overhead
#'   for each additional core. 
#'   
#'   Although this also means that the matrices can be larger than the available
#'   RAM, in practice we find that this slows down the procedure by several 
#'   orders of magnitude. For optimal performance, there should be sufficient
#'   memory to load in each gene expression, correlation, and network matrix
#'   for each dataset. Note: most of this memory is cached; matrices are only 
#'   loaded into RAM when needed (i.e. for the dataset pair for a comparison),
#'   so the physical amount of RAM used will be much lower.
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
#'  of each statistic under the null hypothesis. Two models for the null 
#'  hypothesis are available. Under "overlap", the default, random sampling is
#'  performed only for the set of variables present in both the \emph{discovery} and
#'  \emph{test} datasets. Alternatively, the argument \code{null} can be set to
#'  "all", in which case random sampling is performed on all variables present in
#'  the \emph{test} dataset.
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
#'      Ritchie, S.C., et al., \emph{A scalable permutation approach reveals 
#'      replication and preservation patterns of network modules.} Cell Systems.
#'      \emph{in review} (2016).
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
#'  
#'  For example, \code{results[[1]][[2]][["p.values"]]} is the matrix of 
#'  module preservation p-values when assessing the preservation of modules from
#'  dataset 1 in dataset 2. If \code{simplify = TRUE} then the list structure 
#'  will be simplified where possible.
#'  
#' @seealso 
#'   Functions for: 
#'   \link[=bigMatrix]{bigMatrix objects},
#'   \link[=plotModule]{visualising network modules},
#'   \link[=networkProperties]{calculating module topology}, 
#'   \link[=permutationTest]{calculating permutation test P-values}, and 
#'   \link[=combineAnalyses]{splitting computation over multiple machines}.
#'   
#' @examples
#' \dontrun{
#' # load in example data, correlation, and network matrices for a discovery and test dataset:
#' data("NetRep")
#' 
#' # Convert them to the 'bigMatrix' format:
#' discovery_data <- as.bigMatrix(discovery_data)
#' discovery_correlation <- as.bigMatrix(discovery_correlation)
#' discovery_network <- as.bigMatrix(discovery_network)
#' test_data <- as.bigMatrix(test_data)
#' test_correlation <- as.bigMatrix(test_correlation)
#' test_network <- as.bigMatrix(test_network)
#' 
#' # Set up input lists for each input matrix type across datasets:
#' data_list <- list(discovery=discovery_data, test=test_data)
#' correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
#' network_list <- list(discovery=discovery_network, test=test_network)
#' labels_list <- list(discovery=module_labels)
#' 
#' # Assess module preservation.
#' preservation <- modulePreservation(
#'  data=data_list, correlation=correlation_list, network=network_list,
#'  moduleAssignments=labels_list, nPerm=1000, discovery="discovery", 
#'  test="test"
#' )
#' 
#' }
#' 
#' @import foreach
#' @export
modulePreservation <- function(
  data=NULL, correlation, network, moduleAssignments, modules=NULL, 
  backgroundLabel="0", discovery=1, test=2, selfPreservation=FALSE,
  nCores=NULL, nPerm=NULL, null="overlap", alternative="greater", 
  statCorMethod="pearson", simplify=TRUE, verbose=TRUE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  # Identify the null hypothesis to use 
  nullModels <- c("overlap", "all")
  if (is.na(pmatch(null, nullModels))) {
    stop("overlap must match one of ", paste(nullModels, collapse=" "))
  }
  model <- pmatch(null, nullModels)
  
  # Identify the alternate hypothesis to use
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  alternative <- validAlts[altMatch]
  
  # Set up the correlation measure to use
  validMethods <- c("pearson", "spearman", "kendall", "bicor")
  if (statCorMethod %nin% validMethods)
    stop("'statCorMethod' must be one of ", paste(validMethods, collpase=" "))
  
  # Check for WGCNA if method == "bicor"
  if (statCorMethod == "bicor") {
    tryCatch({
      sink(file.path(tmp.dir, "suppressedWGCNAstartupMessage.txt"))
      suppressMessages(suppressWarnings(requireNamespace("WGCNA")))
      cor <- function(...) WGCNA::bicor(..., quick=1, nThreads=1)[,]
    }, error = function(e) {
      stop("'statCorMethod' = 'bicor' requires the 'WGCNA' package to be",
           "installed.")
    }, finally = { 
      sink() 
    })
  } else {
    cor <- function(...) stats::cor(..., method=statCorMethod)[,]
  }
  
  # Validate 'nPerm'. If 'NULL', we need to process the rest of the input to
  # determine.
  if (!is.null(nPerm) && (!is.numeric(nPerm) || length(nPerm) > 1 || nPerm < 1)) {
    stop("'nPerm' must be a single number > 1")
  }
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=TRUE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # For the other functions this isn't necessary.
  if (missing(moduleAssignments))
    moduleAssignments <- NULL
  if (is.null(moduleAssignments))
    stop("'moduleAssignments' must be provided for each 'discovery' dataset")
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  discovery <- finput$discovery
  test <- finput$test
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)

  # If NULL, automatically determine.
  if (is.null(nPerm)) {
    # If missing, set as the required number for Bonferroni correction.
    # Bonferonni correct for the total number of modules, multiplied the number
    # of datasets each module is tested in.
    multiplier <- sum(sapply(modules, length) * sapply(test, length))
    nPerm <- max(1000, requiredPerms(0.05/multiplier))
  }

  vCat(verbose, 0, "User input ok!")
  
  # Set up return list 
  res <- rep(list(NULL), nDatasets)
  names(res) <- datasetNames
  res <- lapply(res, function(x) { 
    l <- rep(list(NULL), nDatasets)
    names(l) <- datasetNames
    l
  })
  
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
          " from dataset ", '"', di, '"', " within dataset ", '"', di, '"', "."
        )
        next
      }
      tryCatch({
        vCat(
          verbose, 0, sep="", 
          "Calculating preservation of network subsets from dataset ", '"', 
          di, '"', " in dataset ", '"', ti, '"', "."
        )

        # Calculate the overlap between datasets
        ct <- contingencyTable(moduleAssignments, modules, network, di, ti)
        contingency <- ct$contingency
        propVarsPres <- ct$propVarsPres
        overlapVars <- ct$overlapVars
        varsPres <- ct$varsPres
        overlapModules <- ct$overlapModules
        overlapAssignments <- ct$overlapAssignments
        
        nStatistics <- ifelse(!is.null(scaledData[[di]]), 7, 4)
        nModules <- length(overlapModules)
        
        #---------------------------------------------------------------------
        # Get the network properties for each module in the test dataset.
        #   - We only need to calculate this once
        #---------------------------------------------------------------------
        vCat(verbose, 1, "Calculating observed test statistics...")
        discProps <- rep(list(NULL), nModules)
        names(discProps) <- overlapModules
        for (mi in overlapModules) {
          modVars <- names(overlapAssignments %sub_in% mi)
          modInds <- match(modVars, colnames(correlation[[di]]))
          discProps[[mi]] <- moduleProps(network[[di]], modInds, scaledData[[di]])
        }
        
        #---------------------------------------------------------------------
        # Calculate the observed value for each test statistic
        #---------------------------------------------------------------------
        observed <- matrix(NA, nrow=nModules, ncol=nStatistics)
        rownames(observed) <- overlapModules
        for (mi in overlapModules) {
          modVars <- names(overlapAssignments %sub_in% mi)
          discInds <- match(modVars, colnames(correlation[[di]]))
          testInds <- match(modVars, colnames(correlation[[ti]]))
          testProps <-  moduleProps(network[[ti]], testInds, scaledData[[ti]])
          stats <- calcStats(
            discProps[[mi]], testProps, 
            correlation[[di]], discInds,
            correlation[[ti]], testInds
          )
          observed[mi,] <- stats
        }
        colnames(observed) <- names(stats)
        
        # Clean up memory 
        rm(testProps, stats)
        gc()
        
        #---------------------------------------------------------------------
        # Run permutation procedure
        #---------------------------------------------------------------------
        vCat(
          verbose, 1, "Calculating null distributions with", nPerm, 
          "permutations..."
        )
        if(verbose) {
          # To log progress, we will write our progress to a file for each chunk
          while (TRUE) {
            run.dir <- file.path(tmp.dir, paste0(".run-progress", getUUID()))
            # Handle the infintesimally small chance of a UUID collision
            tryCatch({
              dir.create(run.dir)
              break
            }, warning = function(w) {
              if(!grepl("already exists", w$message)) {
                break
              }
            })
          }  
          on.exit({
            unlink(run.dir, recursive=TRUE)
          }, add=TRUE)
        }
        
        foreach(chunk=ichunkTasks(verbose, nPerm, nCores)) %dopar% {
          if (verbose && length(chunk) == 1 && chunk == -1) {
            monitorProgress(nCores - 1, 2, run.dir)
            NULL
          } else {
            tryCatch({
              if (verbose) {
                conns <- setupParProgressLogs(chunk, nCores - 1, 2, run.dir)
                progressBar <- conns[[1]]
              } 
              
              # Attach matrices.
              if (!is.null(scaledData[[di]]))
                scaledData[[di]] <- attach.bigMatrix(scaledData[[di]])
              correlation[[di]] <- attach.bigMatrix(correlation[[di]])
              network[[di]] <- attach.bigMatrix(network[[di]])
              
              if (!is.null(scaledData[[ti]]))
                scaledData[[ti]] <- attach.bigMatrix(scaledData[[ti]])
              correlation[[ti]] <- attach.bigMatrix(correlation[[ti]])
              network[[ti]] <- attach.bigMatrix(network[[ti]])
              
              #---------------------------------------------------------------
              # Calculate the module preservation statistics for each module
              # on a random subset of the same size in the test dataset
              #---------------------------------------------------------------
              chunkStats <- array(
                NA, dim=c(nModules, nStatistics, length(chunk))
              )
              dimnames(chunkStats)[1:2] <- dimnames(observed)
              dimnames(chunkStats)[[3]] <- paste0("permutation.", chunk)
              for (pi in seq_along(chunk)) {
                for (mi in overlapModules) {
                  modVars <- names(overlapAssignments %sub_in% mi)
                  discInds <- match(modVars, colnames(correlation[[di]]))
                  
                  # Select a random subset of nodes of the same size as the 
                  # subsetss, depending on our null model.
                  modSize <- length(modVars)
                  if (model == "overlap") {
                    permVars <- sample(names(overlapVars), modSize)
                  } else {
                    permVars <- sample(colnames(correlation[[ti]]), modSize)
                  }
                  permInds <- match(permVars, colnames(correlation[[ti]]))
                  # Ensure crashes aren't fatal
                  tryCatch({
                    permProps <- moduleProps(network[[ti]], permInds, scaledData[[ti]])
                    chunkStats[mi,,pi] <- calcStats(
                      discProps[[mi]], permProps, 
                      correlation[[di]], discInds,
                      correlation[[ti]], permInds
                    )
                    rm(permProps)
                    gc()
                  }, error = function(e) {
                    warning(
                      "Calculation for module ", mi, " failed on ",
                      "permutation ", chunk[pi], " with error message:\n",
                      e$message
                    )
                  })
                }
                # Update the progress at the end of the loop.
                if (verbose) {
                  updateParProgress(progressBar, chunk[pi])
                  if (nCores == 1) {
                    reportProgress(2, run.dir)
                    if (chunk[pi] == nPerm) {
                      cat("\n")
                    }
                  }
                }
              }
              chunkNum <- ceiling(chunk[1]/length(chunk))
              permFile <- paste0("chunk", chunkNum, "permutations.rds")
              saveRDS(chunkStats, file.path(tmp.dir, permFile))
            }, finally = {
              if (verbose) {
                lapply(conns, close)
              }
              # detach once finished.
              if (!is.null(scaledData[[di]]))
                scaledData[[di]] <- detach.bigMatrix(scaledData[[di]])
              correlation[[di]] <- detach.bigMatrix(correlation[[di]])
              network[[di]] <- detach.bigMatrix(network[[di]])
              
              if (!is.null(scaledData[[ti]]))
                scaledData[[ti]] <- detach.bigMatrix(scaledData[[ti]])
              correlation[[ti]] <- detach.bigMatrix(correlation[[ti]])
              network[[ti]] <- detach.bigMatrix(network[[ti]])
              gc()
            })
          }
        }
        #---------------------------------------------------------------------
        # Load in permutation results
        #---------------------------------------------------------------------
        rm(discProps)
        gc()
        
        nulls <- array(NA, dim=c(nModules, nStatistics, nPerm))
        dimnames(nulls)[[3]] <- rep("", dim(nulls)[3])
        chunkFiles <- list.files(tmp.dir, "chunk[0-9]*permutations.rds")
        offset <- 1
        for (cf in chunkFiles) {
          chunk <- readRDS(file.path(tmp.dir, cf))
          nCPerm <- dim(chunk)[3]
          nulls[,,offset:(offset+nCPerm-1)] <- chunk
          dimnames(nulls)[1:2] <- dimnames(chunk)[1:2]
          dimnames(nulls)[[3]][offset:(offset+nCPerm-1)] <- dimnames(chunk)[[3]]
          offset <- offset + nCPerm
        }
        
        #----------------------------------------------------------------------
        # Order statistics: First density stats, then connectivity, then hybrid
        #----------------------------------------------------------------------
        if (!is.null(scaledData[[di]])) {
          statOrder <- c(
            "avg.weight", "coherence", 
            "cor.cor", "cor.degree", "cor.contrib",
            "avg.cor", "avg.contrib"
          ) 
        } else {
          statOrder <- c("avg.weight", "cor.degree", "cor.cor", "avg.cor")
        }
        observed <- observed[, statOrder, drop=FALSE]
        nulls <- nulls[, statOrder, , drop=FALSE]
        
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
          "Failed with error:\n", e$message, "\nSkipping to next comparison"
        )
      })
    }
  }
  
  # Simplify the output data structure where possible
  for (di in rev(seq_along(res))) {
    for (ti in rev(seq_along(res[[di]]))) {
      if (is.null(res[[di]][[ti]])) {
        res[[di]][[ti]] <- NULL
      }
    }
    if (is.null(res[[di]]) || length(res[[di]]) == 0) {
      res[[di]] <- NULL
    }
  }

  if (simplify) {
    res <- simplifyList(res, depth=2)
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
  res
}
