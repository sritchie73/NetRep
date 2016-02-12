#' Replication and preservation of network modules across datasets
#' 
#' Quantify the preservation of network modules (sub-graphs) in an independent
#' dataset through permutation testing on module topology. Seven network
#' statistics (see details) are calculated for each module and then tested by
#' comparing to distributions generated from their calculation on random subsets
#' in the test dataset.
#'
#' @inheritParams common_params
#' 
#' @param discovery name or index denoting the discovery dataset(s).
#' @param test name or index denoting the dataset(s) to test module preservation 
#'  in.
#' @param nPerm number of permutations to use. Can be specified as a vector if
#'  a different number of permutations is desired for each discovery dataset.
#'  If not specified, the number of permutations will be automatically 
#'  determined (see details).
#' @param nCores number of cores to parallelise the permutation procedure over.
#'  Ignored if the user has already registered a parallel backend.
#' @param exclude an optional vector of modules to exclude from the analysis. If
#'   there are multiple discovery datasets a list of vectors may be provided.
#' @param include an optional vector of modules to include in the
#'   analysis. If there are multiple discovery datasets a list of vectors may be
#'   provided.
#' @param null variables to include when generating the null distributions. 
#'  Must be either "overlap" or "all" (see details).
#' @param alternative The type of module preservation test to perform. Must be 
#'   one of "greater" (default), "less" or "two.sided" (see details).
#' @param statCorMethod character vector indicating method to use when calculating 
#'   the correlation based statistics (see details). Must be one of "pearson", 
#'   "spearman", or "kendall". If the WGCNA package is installed then "bicor" 
#'   may also be specified as an option (see \code{\link[WGCNA]{bicor}}). 
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#' @param simplify logical; if \code{TRUE}, simplify the structure of the output
#'  list if possible (see Return Value).
#' @param keepNulls logical; if \code{TRUE}, the null distributions are returned
#'  as part of the output.
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
#'   \link[=bigMatrix-get]{converting to, loading in}, and 
#'   \link[=bigMatrix-out]{writing out} 'bigMatrix' objects.
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
#'  calculated from the module’s \emph{summary profile}. This is the eigenvector
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
#'  \code{\link{plotNodeContribution}}.) Similarly, a strong
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
#'  non-meaningful, variations in weighted degree.
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
#'  be accurate (see \code{\link{perm.test}}) \emph{(3)}.
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
#'  The returned data structure is organised as a nested list of lists, which 
#'  should be accessed as \code{results[[discovery]][[test]]}. If
#'  \code{simplify} is set to \code{TRUE}, then this structure will be 
#'  simplified as much as possible depending on the combination of dataset 
#'  comparisons that have been performed. For each dataset-comparison a list of
#'  the following objects are returned:
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
#'      and \emph{test} datasets, then a contigency table showing the overlap
#'      between modules across datasets is returned. Rows correspond to modules
#'      in the \emph{discovery} dataset, columns to modules in the \emph{test}
#'      dataset.
#'    }
#'  }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' 
#' ## Example 1: Assess replication of all modules from one 
#' ## cohort in an independent dataset
#' 
#' # First we need some example data
#' geA <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geA) <- paste0("Gene_", 1:100)
#' rownames(geA) <- paste0("CohortA_", 1:50)
#' coexpA <- cor(geA) # correlation
#' adjA <- abs(coexpA)^5 # adjacency
#' moduleAssignments <- sample(1:7, size=100, replace=TRUE)
#' names(moduleAssignments) <- paste0("Gene_", 1:100)
#' 
#' geB <- matrix(rnorm(70*100), ncol=100) # gene expression
#' colnames(geB) <- paste0("Gene_", 1:100) 
#' rownames(geB) <- paste0("CohortB_", 1:70)
#' coexpB <- cor(geB) # correlation
#' adjB <- abs(coexpB)^6 # adjacency
#' 
#' # Now format the data for input to modulePreservation
#' data <- list(
#'   cohortA=as.bigMatrix(geA, "geA_bm"),
#'   cohortB=as.bigMatrix(geA, "geB_bm")    
#' )
#' correlation <- list(
#'   cohortA=as.bigMatrix(coexpA, "coexpA_bm"),
#'   cohortB=as.bigMatrix(coexpB, "coexpB_bm")
#' )
#' adjacency <- list(
#'   cohortA=as.bigMatrix(adjA, "adjA_bm"),
#'   cohortB=as.bigMatrix(adjB, "adjB_bm")
#' )
#' 
#' # Assess module preservation, using two cores
#' replication <- modulePreservation(
#'   data, correlation, adjacency, moduleAssignments, 
#'   nCores=2
#' )
#' 
#' ## Example 2: assess replication of two disease-associated modules
#' replication <- modulePreservation(
#'   data, correlation, adjacency, moduleAssignments,
#'   nCores=2, include = c("4", "7")
#' )
#' 
#' ## Example 3: exclude a module from the analysis
#' replication <- modulePreservation(
#'   data, correlation, adjacency, moduleAssignments,
#'   nCores=2, exclude = "0"
#' )
#' 
#' ## Example 4: assess preservation of modules across multiple
#' ## tissues
#' geAdipose <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geAdipose) <- paste0("Gene_", 1:100)
#' rownames(geAdipose) <- paste0("Sample_", 1:50)
#' coexpAdipose <- cor(geAdipose) # correlation
#' adjAdipose <- abs(coexpAdipose)^5 # adjacency
#' adiposeModules <- sample(0:7, size=100, replace=TRUE)
#' names(adiposeModules) <- paste0("Gene_", 1:100)
#' 
#' geLiver <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geLiver) <- paste0("Gene_", 1:100)
#' rownames(geLiver) <- paste0("Sample_", 1:50)
#' coexpLiver <- cor(geLiver) # correlation
#' adjLiver <- abs(coexpLiver)^6 # adjacency
#' liverModules <- sample(0:12, size=100, replace=TRUE)
#' names(liverModules) <- paste0("Gene_", 1:100)
#' 
#' geHeart <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geHeart) <- paste0("Gene_", 1:100)
#' rownames(geHeart) <- paste0("Sample_", 1:50)
#' coexpHeart <- cor(geHeart) # correlation
#' adjHeart <- abs(coexpHeart)^4 # adjacency
#' heartModules <- sample(0:5, size=100, replace=TRUE)
#' names(heartModules) <- paste0("Gene_", 1:100)
#' 
#' # Now format the data for input to modulePreservation
#' data <- list(
#'   adipose=as.bigMatrix(geAdipose, "geAdipose_bm"),
#'   liver=as.bigMatrix(geLiver, "geLiver_bm"),  
#'   heart=as.bigMatrix(geHeart, "geHeart_bm") 
#' )
#' correlation <- list(
#'   adipose=as.bigMatrix(coexpAdipose, "coexpAdipose_bm"),
#'   liver=as.bigMatrix(coexpLiver, "coexpLiver_bm"),  
#'   heart=as.bigMatrix(coexpHeart, "coexpHeart_bm") 
#' )
#' adjacency <- list(
#'   adipose=as.bigMatrix(adjAdipose, "adjAdipose_bm"),
#'   liver=as.bigMatrix(adjLiver, "adjLiver_bm"),  
#'   heart=as.bigMatrix(adjHeart, "adjHeart_bm") 
#' )
#' moduleAssignments <- list(
#'   adipose=adiposeModules, liver=liverModules, heart=heartModules
#' )
#' 
#' # Assess the preservation of each module in each non-discovery
#' # tissue.
#' preservation <- modulePreservation(
#'   data, correlation, adjacency, moduleAssignments,
#'   nCores=2, discovery=c("adipose", "liver", "heart"), 
#'   test=c("adipose", "liver", "heart")
#' )
#' 
#' # Remove bigMatrix files used in examples
#' unlink("*_bm*")
#' }
#' 
#' @import foreach
#' @import RhpcBLASctl
#' @export
modulePreservation <- function(
  data=NULL, correlation, network, moduleAssignments,
  discovery=1, test=2, nCores=1, nPerm, exclude,
  include, null="overlap", alternative="greater", statCorMethod="pearson",
  simplify=TRUE, verbose=TRUE, keepNulls=FALSE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  # Temporary directory to store new bigMatrix objects in
  vCat(verbose, 0, "creating directory for temporary objects...")
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    vCat(verbose, 0, "removing temporary object directory...")
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
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
  
  if (!is.numeric(nCores) || length(nCores) > 1 || nCores < 1)
    stop("'nCores' must be a single number greater than 0")
  if (!is.null(data) && !is.list(data))
    stop("Expecting a list of matrices for argument 'data'")
  if (!is.list(correlation))
    stop("Expecting a list of matrices for argument 'correlation'")
  if (!is.list(network))
    stop("Expecting a list of matrices for argument 'network'")
  
  # Check for valid corMethod options
  validMethods <- c("pearson", "spearman", "kendall", "bicor")
  if (statCorMethod %nin% validMethods)
    stop("'corMethod' must be one of ", paste(validMethods, collpase=" "))
  
  # Check for WGCNA. We need to open a temporary sink to suppress all of WGCNA's
  # startup messages.
  hasWGCNA <- FALSE
  sink(file.path(tempdir(), "suppressedWGCNAstartupMessage.txt"))
  if (suppressMessages(suppressWarnings(requireNamespace("WGCNA"))))
    hasWGCNA <- TRUE
  sink()
  if (statCorMethod == "bicor" & !hasWGCNA)
    stop("statCorMethod='bicor' requires the WGCNA package.")
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structure.
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(correlation), names(correlation)
  )
  
  # Try to intelligently handle different types of user input
  data <- dynamicMatLoad(data)
  correlation <- dynamicMatLoad(correlation)
  network <- dynamicMatLoad(network)
  
  # Sanity check input for consistency.
  checkSets(
    data, correlation, network, moduleAssignments, discovery, test
  )
  
  # Are the datasets named, or will we iterate by index?
  if (!is.null(names(correlation))) {
    datasets <- names(correlation)
    if (class(discovery) != "character")
      discovery <- datasets[discovery]
    if (class(test) != "character")
      test <- datasets[test]
  } else {
    datasets <- seq_along(correlation)
  }
  nDatasets <- length(datasets)
  
  # same for the include and exclude modules arguments
  include <- formatInclude(
    include, discovery, length(correlation), datasets
  )
  exclude <- formatExclude(
    exclude, discovery, length(correlation), datasets
  )
  
  # Sanity check input data for values that will cause the calculation of 
  # network properties and statistics to hang.
  vCat(verbose, 0, "checking matrices for non-finite values...")
  lapply(data, checkFinite)
  lapply(correlation, checkFinite)
  lapply(network, checkFinite)
  
  # Temporarily create scaled data set for the calculation of the
  # summary expression profile
  sge <- NULL
  if (!is.null(data))
    sge <- lapply(data, scaleBigMatrix, tmp.dir)
  
  # Set up return list 
  res <- rep(list(NULL), nDatasets)
  res <- lapply(res, function(x) { 
    l <- rep(list(NULL), nDatasets)
    if (class(datasets) == "character")
      names(l) <- datasets
    l
  })
  if (class(datasets) == "character")
    names(res) <- datasets
  
  #-----------------------------------------------------------------------------
  # Atuomatically determine the number of permutations if unspecified
  #-----------------------------------------------------------------------------
  if (missing(nPerm)) {
    # If missing, set as the required number for Bonferroni correction.
    # Bonferonni correct for the total number of modules, multiplied the number
    # of datasets each module is tested in.
    multiplier <- sum(sapply(discovery, function(di) {
      modules <- names(table(moduleAssignments[[di]]))
      if (!is.null(exclude[[di]])) {
        modules <- modules %sub_nin% exclude[[di]]
      }
      if (!is.null(include[[di]])) {
        modules <- modules %sub_in% include[[di]]
      }
      nTest <- length(test %sub_nin% di)
      length(modules)*nTest
    }))
    nPerm <- max(1000, requiredPerms(0.05/multiplier))
    nPerm <- rep(nPerm, length(discovery))
    names(nPerm) <- discovery
  } else if (length(nPerm) > 1) {
    if (!is.numeric(nPerm))
      stop("'nPerm' must be a numeric vector")
    if (is.null(names(nPerm))) {
      names(nPerm) <- discovery
    } else if (any(names(nPerm) %nin% discovery)) {
      stop("mismatch between dataset names in 'nPerm' and 'discovery'")
    }
  } else if (!is.numeric(nPerm) ){
    stop("'nPerm' must be numeric")
  } else {
    nPerm <- rep(nPerm, length(discovery))
    names(nPerm) <- discovery
  }
  
  #-----------------------------------------------------------------------------
  # Set up parallel backend
  #-----------------------------------------------------------------------------
  
  # First, check whether the user has already set up a parallel backend. In this
  # case, we can ignore the `nCores` argument.
  if (getDoParWorkers() > 1) {
    vCat(
      verbose, 0, "Using user-registered parallel backend with 1 reporter core",
      "and", getDoParWorkers() - 1, "worker cores."
    )
    nCores <- getDoParWorkers()
  } 
  # If the user is on a Windows machine, we have to use the `doParallel` package 
  else if (.Platform$OS.type == "windows" & nCores > 1) {
    # Quietly load parallel backend packages. Throw our own warning and 
    # continue
    if(suppressWarnings(suppressMessages(requireNamespace("doParallel")))) {
      # we need an additional thread to monitor and report progress
      if (verbose)  
        nCores <- nCores + 1
      cl <- parallel::makeCluster(nCores)
      doParallel::registerDoParallel(cl)
      on.exit({
        parallel::stopCluster(cl)
      }, add=TRUE)
      vCat(verbose, 0, "Running on", nCores - 1, "cores.")
      if ((nCores - 1) > parallel::detectCores()) {
        stop(
          "Requested number of threads (", nCores - 1, ") is higher than the ",
          "number of available cores (", parallel::detectCores(), 
          "). Using too many threads may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end 
      # once the analysis has finished.
      vCat(
        TRUE, 0, file=stderr(),
        "Warning: unable to find 'doParallel' package, running on 1 core." 
      )
    }
  } else if (.Platform$OS.type == "unix" & nCores > 1) {
    # Quietly load parallel backend packages. Throw our own warning and 
    # continue
    if(suppressWarnings(suppressMessages(requireNamespace("doMC")))) {
      # we need an additional thread to monitor and report progress
      if (verbose) 
        nCores <- nCores + 1
      doMC::registerDoMC(nCores)
      vCat(verbose, 0, "Running on", nCores - 1, "cores.")
      if ((nCores - 1) > parallel::detectCores()) {
        stop(
          "Requested number of threads (", nCores - 1, ") is higher than the ",
          "number of available cores (", parallel::detectCores(), 
          "). Using too many threads may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end 
      # once the analysis has finished.
      vCat(
        TRUE, 0, file=stderr(),
        "Unable to find 'doMC' package, running on 1 core."
      )
    }
  } else {
    vCat(verbose, 0, "Running on 1 cores.")
  }
  
  # Suppress annoying foreach warning generated when using %dopar% and running 
  # in serial
  if (nCores == 1) {
    suppressWarnings({
      ii <- 0 # suppress R CMD check note
      foreach(ii = 1:2) %dopar% { ii }
    })
  }
  
  # Since we expect the user to explicitly handle the number of parallel threads,
  # we will disable the potential implicit parallelism on systems where R has
  # been compiled against a multithreaded BLAS, e.g. OpenBLAS. 
  omp_set_num_threads(1)
  blas_set_num_threads(1)
  
  #-----------------------------------------------------------------------------
  # Set up correlation function
  #-----------------------------------------------------------------------------
  if (statCorMethod == "bicor") {
    cor <- function(...) WGCNA::bicor(..., quick=1, nThreads=1)[,]
  } else {
    cor <- function(...) stats::cor(..., method=statCorMethod)[,]
  }
  
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
  for (di in datasets) {
    for (ti in datasets) {
      if ((di %in% discovery) & (ti %in% test) & (di != ti)) {
        tryCatch({
          vCat(
            verbose, 0, sep="", 
            "Calculating preservation of network subsets from dataset ",
            di, " in dataset ", ti, "."
          )
          #---------------------------------------------------------------------
          # Set up variables for this comparison
          #---------------------------------------------------------------------
          
          # Force modules to be character vectors to avoid improper indexing
          if (is.numeric(moduleAssignments[[di]])) { 
            moduleAssignments[[di]] <- as.character(moduleAssignments[[di]])
            names(moduleAssignments[[di]]) <- colnames(correlation[[di]])
          }
          
          # To simplify later function calls, we need to get a vector of module
          # assignments only for (a) modules of interest and (b) the variables
          # present in both datasets for those modules.
          overlapVars <- intersect(
            colnames(correlation[[di]]), 
            colnames(correlation[[ti]])
          )
          overlapAssignments <- moduleAssignments[[di]][overlapVars]
          # Restrict to modules of interest
          modules <- unique(moduleAssignments[[di]])
          if (!is.null(exclude[[di]]))
            modules <- modules %sub_nin% exclude[[di]]
          if (!is.null(include[[di]]))
            modules <- modules %sub_in% include[[di]]
          overlapAssignments <- overlapAssignments %sub_in% modules 
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
          modulesWithNoOverlap <- modules %sub_nin% overlapAssignments
          varsPres <- c(varsPres, rep(0, length(modulesWithNoOverlap)))
          names(varsPres)[names(varsPres) == ""] <- modulesWithNoOverlap
          varsPres <- varsPres[orderAsNumeric(names(varsPres))]
          
          # What proportion?
          moduleSizes <- table(moduleAssignments[[di]])
          moduleSizes <- moduleSizes[names(moduleSizes) %sub_in% modules]
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
            contingency <- contingency[modules,,drop=FALSE]
            
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
          
          nStatistics <- ifelse(!is.null(sge[[di]]), 7, 4)
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
            discProps[[mi]] <- moduleProps(network[[di]], modInds, sge[[di]])
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
            testProps <-  moduleProps(network[[ti]], testInds, sge[[ti]])
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
            verbose, 1, "Calculating null distributions with", nPerm[di], 
            "permutations..."
          )
          if(verbose) {
            # To log progress, we will write our progress to a file for each chunk
            while (TRUE) {
              run.dir <- file.path(tempdir(), paste0(".run-progress", getUUID()))
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
          foreach(chunk=ichunkTasks(verbose, nPerm[di], nCores)) %dopar% {
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
                if (!is.null(data))
                  sge <- lapply(sge, attach.bigMatrix)
                correlation <- lapply(correlation, attach.bigMatrix)
                network <- lapply(network, attach.bigMatrix)
                
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
                    
                    # Select a random subset of nodes of the same size as the subset 
                    # ss, depending on our null model.
                    modSize <- length(modVars)
                    if (model == "overlap") {
                      permVars <- sample(names(overlapVars), modSize)
                    } else {
                      permVars <- sample(colnames(correlation[[ti]]), modSize)
                    }
                    permInds <- match(permVars, colnames(correlation[[ti]]))
                    # Ensure crashes aren't fatal
                    tryCatch({
                      permProps <- moduleProps(network[[ti]], permInds, sge[[ti]])
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
                      if (chunk[pi] == nPerm[di]) {
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
                if (!is.null(data))
                  sge <- lapply(sge, detach.bigMatrix)
                correlation <- lapply(correlation, detach.bigMatrix)
                network <- lapply(network, detach.bigMatrix)
              })
            }
          }
          #---------------------------------------------------------------------
          # Load in permutation results
          #---------------------------------------------------------------------
          rm(discProps)
          gc()
          
          nulls <- array(NA, dim=c(nModules, nStatistics, nPerm[di]))
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
          
          #---------------------------------------------------------------------
          # Calculate permutation p-value
          #---------------------------------------------------------------------
          vCat(verbose, 1, "Calculating P-values...")
          p.values <- matrix(NA, nrow=nModules, ncol=nStatistics)
          dimnames(p.values) <- dimnames(observed)
          for (mi in overlapModules) {
            for (si in seq_len(nStatistics)) {
              # Does the order of nodes in each permutation affect the statistic?
              if (colnames(observed)[si] %in% c("avg.weight", "coherence")) {
                order <- FALSE
              } else {
                order <- TRUE
              }
              
              # Get the p-values
              p.values[mi, si] <- perm.test(
                nulls[mi, si, ], observed[mi, si], 
                varsPres[mi], length(overlapVars),
                order=order, alternative=alternative
              )
              
            }
          }
          
          # Order statistics: First density stats, then connectivity, then hybrid
          if (!is.null(sge[[di]])) {
            statOrder <- c(
              "avg.weight", "coherence", 
              "cor.cor", "cor.degree", "cor.contrib",
              "avg.cor", "avg.contrib"
            ) 
          } else {
            statOrder <- c("avg.weight", "cor.degree", "cor.cor", "avg.cor")
          }
          
          
          #---------------------------------------------------------------------
          # Collate results
          #---------------------------------------------------------------------
          vCat(verbose, 1, "Collating results...")
          if (!keepNulls) {
            nulls <- NULL # ha!
          } else {
            nulls <- nulls[, statOrder,]
          }
          res[[di]][[ti]] <- list(
            observed = observed[, statOrder],
            nulls = nulls,
            p.values = p.values[, statOrder],
            nVarsPresent = varsPres,
            propVarsPresent = propVarsPres,
            contingency = contingency
          )
          # remove NULL outputs
          res[[di]][[ti]] <- res[[di]][[ti]][
            !sapply(res[[di]][[ti]], is.null)
          ]
          
          gc()
          vCat(verbose, 0, "Done!")
        }, error=function(e) {
          warning(
            "Failed with error:\n", e$message, "\nSkipping to next comparison"
          )
        })
      }
    }
  }
  # Simplify the output data structure where possible
  if (simplify) {
    if (length(discovery) == 1 && length(test) == 1) {
      res <- res[[discovery]][[test]]
    } else if (length(discovery) == 1) {
      res <- res[[discovery]][test] 
    } else if (length(test) == 1) {
      res <- lapply(res[discovery], `[[`, test)
    }
  }
  res
}
