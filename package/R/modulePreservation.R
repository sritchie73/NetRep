#' Replication and preservation of network modules across datasets
#' 
#' Assess whether gene coexpression modules replicate or are preserved in an 
#' independent dataset. Module preservation is assessed using a permutation 
#' procedure performed on seven module preservation statistics (see details). 
#' Each dataset requires a precalculated coexpression network and adjecency 
#' network (see details). Providing the gene expression is optional, but 
#' recommended. Matrix data is ideally stored in the
#' \code{\link[=bigMatrix-class]{bigMatrix}} format (see details).
#' 
#' @param geneExpression A list of gene expression matrices, one for each 
#' dataset. Expects columns to be genes, rows to be samples. 
#' @param coexpression A list of coexpression matrices, one for each dataset.
#' @param adjacency A list of adjacency matrices, one for each dataset.
#' @param moduleAssignments A list of vectors for each \emph{discovery} dataset
#'  containing the module assignments for each gene in the respective dataset.
#' @param discovery datasets where module discovery has been performed.
#' @param test datasets to test for preservation of modules from each 
#'  \emph{discovery} dataset.
#' @param nPerm number of permutations to use. Can be specified as a vector if
#'  a different number of permutations is desired for each discovery dataset.
#'  If not specified, the number of permutations will be automatically 
#'  determined (see details).
#' @param nCores number of cores to parallelise the permutation procedure over.
#' @param lowmem logical; should memory usage be minimised? Useful on machines
#'  with limited RAM when running in parallel.
#' @param excludeModules optional list of vectors containing modules to exclude 
#'  from the analysis for each \code{discovery} dataset. If unspecified, the 
#'  preservation of all modules will be tested.
#' @param includeModules optional list of vectors containing modules to include 
#'  in the analysis for each \code{discovery} dataset. If unspecified, the 
#'  preservation of all modules will be tested.
#' @param null the type of null model, either "overlap" or "all" (see details).
#' @param alternative The type of module preservation test to perform. Must be 
#'   one of "greater" (default), "less" or "two.sided" (see details).
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#' @param simplify logical; if \code{TRUE}, simplify the structure of the output
#'  list if possible (see Return Value).
#' @param keepNulls logical; if \code{TRUE}, the null distributions are returned
#'  as part of the output.
#'  
#' @details
#'  \subsection{Input data structure:}{
#'   The topological properties used to assess module preservation are designed 
#'   for networks constructed using Weighted Gene Coexpression Network Analysis 
#'   (\pkg{\link[WGCNA]{WGCNA}}, \emph{(3)}). These are calculated from the gene expression
#'   for each dataset, the pairwise correlation between genes (coexpression) for
#'   each dataset, and the pairwise gene adjacencies (adjacency) for each
#'   dataset. The adjacency is typically the absolute value of the correlation
#'   raised to a power to penalise weak correlations \emph{(3)}. Module
#'   preservation can also be assessed on networks without the gene expression
#'   data, but only a limited subset of the statistics will be calculated.
#'   Network modules are usually clusters of tightly coexpressed genes
#'   \emph{(3)}, but the procedure is also useful for assessing known gene sets,
#'   i.e. pathways across conditions or tissues \emph{(1)}.
#'   
#'   The arguments \code{geneExpression}, \code{coexpression}, and 
#'   \code{adjacency} each expect a \code{\link{list}} where each element 
#'   contains the matrix data for each respective dataset. This matrix data 
#'   should be stored as a 'bigMatrix' object (see 
#'   \link[=bigMatrix-get]{converting matrix data to 'bigMatrix' data}).
#'   Similarly, the \code{moduleAssignments} argument expects a list of named
#'   vectors, which contain the the module assignments for each gene in the
#'   respective dataset. List elements corresponding to datasets where module
#'   discovery has not been performed should contain \code{NULL}, unless the
#'   datasets are named throughout the function arguments. I.e. where the
#'   \code{\link{names}} of \code{geneExpression}, \code{coexpression}, and
#'   \code{adjacency} correspond to the names of each dataset of interest, the
#'   names of the \code{discovery} dataset can be used to look up the respective
#'   module assignments in the \code{moduleAssignments} list. If module
#'   discovery has only been performed in one dataset, then the
#'   \code{moduleAssignments} will also accept a named vector.
#'   
#'   The \code{discovery} arguments specifies which dataset the \code{modules} 
#'   of interest were discovered in, and the \code{test} argument specifies 
#'   which dataset to calculate the network properties from. These arguments are
#'   ignored if data is provided for only one dataset. Otherwise, the function
#'   defaults to calculating the network properties for \code{modules} from the 
#'   first dataset specified in the list structure of \code{geneExpression}, 
#'   \code{coexpression}, and \code{adjacency}, in that same dataset.
#' }
#' \subsection{'bigMatrix' vs. 'matrix' input data:}{
#'   Although the function expects \code{\link[=bigMatrix-class]{bigMatrix}}
#'   data, regular 'matrix' objects are also accepted. In this case, the
#'   'matrix' data is temporarily converted to 'bigMatrix' by the function. This
#'   conversion process involves writing out each matrix as a binary file on
#'   disk, which can take a long time for large datasets. It is strongly
#'   recommended for the user to store their data as 'bigMatrix' objects, as the
#'   \link{modulePreservation} function, \link[=plotModule]{plotting} 
#'   \link[=plotTopology]{functions}, \link[=geneOrder]{gene} and 
#'   \link[=sampleOrder]{sample} ordering also expect 'bigMatrix' objects. 
#'   Further, 'bigMatrix' objects have a number of benefits, including 
#'   instantaneous load time from any future R session, and parallel access from
#'   mutliple independent R sessions. Methods are provided for 
#'   \link[=bigMatrix-get]{converting to, loading in}, and 
#'   \link[=bigMatrix-out]{writing out} 'bigMatrix' objects.
#' }
#' \subsection{Memory usage:}{
#'   For most machines there is a trade-off between memory usage per core and
#'   the number of cores that can be utilised in parallel. The \code{lowmem}
#'   argument controls the RAM usage for per core. When \code{TRUE}, the
#'   default, the permutation procedure will unload the RAM after each
#'   permutation, allowing the permutation procedure to be parallelised over
#'   more cores. \code{lowmem} should not be set to \code{FALSE} unless there is
#'   sufficient RAM to load the all matrices into RAM for each dataset
#'   comparison on each core.
#' }
#' \subsection{Module Preservation Statistics:}{
#'  Module preservation is assessed through seven statistics \emph{(1)}:
#'  \enumerate{
#'    \item{\code{mean.adj}:}{
#'      The mean adjacency, or module density, measures how densely connected a 
#'      module is in the \emph{test} dataset.
#'    }
#'    \item{\code{pve}:}{
#'      Short for "the proportion of variance explained in the underlying gene
#'      expression data for the module by its summary expression profile in the
#'      \emph{test} dataset". The summary expression profile is calculated as
#'      the first eigenvector from a principal component analysis on the
#'      module's (scaled) gene expression data. The summary expression profile
#'      is commonly referred to as the "module eigengene (ME)", and abbreviated
#'      as \code{propVarExpl} in \emph{(1)}.
#'    }
#'    \item{\code{cor.coexp}:}{
#'      The correlation of coexpression patterns for a module across the
#'      \emph{discovery} and \emph{test} datasets. It is also referred to as
#'      the "correlation of correlation" and abbreviated as \code{cor.cor} in 
#'      \emph{(1)}.
#'    }
#'    \item{\code{cor.kIM}:}{
#'      The correlation of intramodular connectivity across the \emph{discovery}
#'      and \emph{test} datasets. Intramodular connectivity is quantified as the
#'      sum of adjacency for a gene to all other genes in the module.
#'    }
#'    \item{\code{cor.MM}:}{
#'      The correlation of intramodular module membership across the
#'      \emph{discovery} and \emph{test} datasets. Module membership denotes the
#'      correlation between each gene and the summary expression profile for
#'      that module in the given dataset. The statistic is abbreviated as
#'      \code{cor.kME} in \emph{(1)}.
#'    }
#'    \item{\code{mean.coexp}:}{
#'      The mean sign-aware coexpression measures the average coexpression for a
#'      module in the \emph{test} dataset, multiplied by the sign of the 
#'      coexpression in the \emph{discovery} dataset. It assesses how strongly 
#'      correlated the genes are in the test dataset (in either direction), 
#'      penalising the module if any gene pairs whose correlation flips between
#'      datasets. It is also referred to as the "mean sign-aware correlation"
#'      and abbreviated as \code{mean.cor} in \emph{(1)}.
#'    }
#'    \item{\code{mean.MM}:}{
#'      The mean sign-aware module membership measures the average module 
#'      membership in the \emph{test} dataset, multiplied by the sign of the 
#'      module membership in the \emph{discovery} dataset. It measures how 
#'      coherent the gene expression is in the \emph{test} dataset, penalising 
#'      the module if any genes are differentially expressed compared to the 
#'      module in one, but not both \emph{discovery} and \emph{test} datasets. 
#'      It is also abbreivated as \code{mean.kME} in \emph{(1)}.
#'    }
#'  }
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
#'  module in the \emph{discovery} dataset and \code{nPerm} random gene sets of
#'  the same size in the \emph{test} dataset in order to assess the distribution
#'  of each statistic under the null hypothesis. Two models for the null 
#'  hypothesis are available. Under "overlap", the default, random sampling is
#'  performed only for the set of genes present in both the \emph{discovery} and
#'  \emph{test} datasets. Alternative, the argument \code{null} can be set to
#'  "all", in which case random sampling is performed on all genes present in
#'  the \emph{test} dataset. The latter may be suitable when assessing module
#'  replication across species.
#'   
#'  The number of permutations required for any given significance threshold is 
#'  approximately 1 / the desired significance for one sided tests, and double 
#'  that for two-sided tests. This can be calculated with 
#'  \code{\link{requiredPerms}}. When \code{nPerm} is not specified, the number 
#'  of permutations is automatically calculated as the number required for a 
#'  Bonferroni corrected significance threshold adjusting for the number of 
#'  modules in each \emph{discovery} dataset multiplied by the number of 
#'  \emph{test} datasets. Although assessing the replication of a small number 
#'  of modules calls for very few permutations, we recommend using no fewer than
#'  200 as fewer permutations are unlikely to generate representative null 
#'  distributions. \strong{Note:} the assumptions used by 
#'  \code{\link{requiredPerms}} break down when assessing the preservation of 
#'  very small modules in a very small dataset (e.g. gene sets in a dataset with
#'  less than 100 genes total). However, the reported p-values will still be 
#'  accurate (see \code{\link{perm.test}}) \emph{(2)}.
#' }
#' 
#' @references 
#'   \enumerate{
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
#'    \item{\code{nGenesPresent}:}{
#'      A vector containing the number of genes that are present in the test
#'      dataset for each module.
#'    }
#'    \item{\code{propGenesPresent}:}{
#'      A vector containing the proportion of genes present in the test dataset
#'      for each module. Modules where this is less than 1 should be 
#'      investigated further before making judgements about preservation to 
#'      ensure that the missing genes are not the most connected ones.
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
#' coexpA <- cor(geA) # coexpression
#' adjA <- abs(coexpA)^5 # adjacency
#' moduleAssignments <- sample(1:7, size=100, replace=TRUE)
#' names(moduleAssignments) <- paste0("Gene_", 1:100)
#' 
#' geB <- matrix(rnorm(70*100), ncol=100) # gene expression
#' colnames(geB) <- paste0("Gene_", 1:100) 
#' rownames(geB) <- paste0("CohortB_", 1:70)
#' coexpB <- cor(geB) # coexpression
#' adjB <- abs(coexpB)^6 # adjacency
#' 
#' # Now format the data for input to modulePreservation
#' geneExpression <- list(
#'   cohortA=as.bigMatrix(geA, "geA_bm"),
#'   cohortB=as.bigMatrix(geA, "geA_bm")    
#' )
#' coexpression <- list(
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
#'   geneExpression, coexpression, adjacency, moduleAssignments, 
#'   nCores=2
#' )
#' 
#' ## Example 2: assess replication of two disease-associated modules
#' replication <- modulePreservation(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   nCores=2, includeModules = c("4", "7")
#' )
#' 
#' ## Example 3: exclude a module from the analysis
#' replication <- modulePreservation(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   nCores=2, excludeModules = "0"
#' )
#' 
#' ## Example 4: assess preservation of modules across multiple
#' ## tissues
#' geAdipose <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geAdipose) <- paste0("Gene_", 1:100)
#' rownames(geAdipose) <- paste0("Sample_", 1:50)
#' coexpAdipose <- cor(geAdipose) # coexpression
#' adjAdipose <- abs(coexpAdipose)^5 # adjacency
#' adiposeModules <- sample(0:7, size=100, replace=TRUE)
#' names(adiposeModules) <- paste0("Gene_", 1:100)
#' 
#' geLiver <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geLiver) <- paste0("Gene_", 1:100)
#' rownames(geLiver) <- paste0("Sample_", 1:50)
#' coexpLiver <- cor(geLiver) # coexpression
#' adjLiver <- abs(coexpLiver)^6 # adjacency
#' liverModules <- sample(0:12, size=100, replace=TRUE)
#' names(liverModules) <- paste0("Gene_", 1:100)
#' 
#' geHeart <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geHeart) <- paste0("Gene_", 1:100)
#' rownames(geHeart) <- paste0("Sample_", 1:50)
#' coexpHeart <- cor(geHeart) # coexpression
#' adjHeart <- abs(coexpHeart)^4 # adjacency
#' heartModules <- sample(0:5, size=100, replace=TRUE)
#' names(heartModules) <- paste0("Gene_", 1:100)
#' 
#' # Now format the data for input to modulePreservation
#' geneExpression <- list(
#'   adipose=as.bigMatrix(geAdipose, "geAdipose_bm"),
#'   liver=as.bigMatrix(geLiver, "geLiver_bm"),  
#'   heart=as.bigMatrix(geHeart, "geHeart_bm") 
#' )
#' coexpression <- list(
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
#'   geneExpression, coexpression, adjacency, moduleAssignments,
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
  geneExpression=NULL, coexpression, adjacency, moduleAssignments,
  discovery=1, test=2, nCores=1, lowmem=TRUE, nPerm, excludeModules,
  includeModules, null="overlap", alternative="greater",
  simplify=TRUE, verbose=TRUE, keepNulls=FALSE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  # Temporary directory to store new bigMatrix objects in
  vCat(verbose, 0, "creating directory for temporary objects...")
  tmp.dir <- paste0(".temp-objects", getUUID())
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
  
  # Identiyf the alternate hypothesis to use
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  
  if (!is.numeric(nCores) || length(nCores) > 1 || nCores < 1)
    stop("'nCores' must be a single number greater than 0")
  if (!is.null(geneExpression) && !is.list(geneExpression))
    stop("Expecting a list of gene expression matrices, one for each dataset")
  if (!is.list(coexpression))
    stop("Expecting a list of coexpression matrices, one for each dataset")
  if (!is.list(adjacency))
    stop("Expecting a list of adjacency matrices, one for each dataset")
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structure.
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression)
  )
  
  # Try to intelligently handle different types of user input
  geneExpression <- dynamicMatLoad(geneExpression, backingpath=tmp.dir)
  coexpression <- dynamicMatLoad(coexpression, backingpath=tmp.dir)
  adjacency <- dynamicMatLoad(adjacency, backingpath=tmp.dir)
  
  # Sanity check input for consistency.
  checkSets(
    geneExpression, coexpression, adjacency, moduleAssignments, discovery, test
  )
  
  # Are the datasets named, or will we iterate by index?
  if (!is.null(names(coexpression))) {
    datasets <- names(coexpression)
    if (class(discovery) != "character")
      discovery <- datasets[discovery]
    if (class(test) != "character")
      test <- datasets[test]
  } else {
    datasets <- seq_along(coexpression)
  }
  nDatasets <- length(datasets)
  
  # same for the include and exclude modules arguments
  includeModules <- formatInclude(
    includeModules, discovery, length(coexpression), datasets
  )
  excludeModules <- formatExclude(
    excludeModules, discovery, length(coexpression), datasets
  )
  
  # Sanity check input data for values that will cause the calculation of 
  # network properties and statistics to hang.
  vCat(verbose, 0, "checking matrices for non-finite values...")
  lapply(geneExpression, checkFinite)
  lapply(coexpression, checkFinite)
  lapply(adjacency, checkFinite)
  
  # Temporarily create scaled gene expression set for the calculation of the
  # summary expression profile
  sge <- NULL
  if (!is.null(geneExpression))
    sge <- lapply(geneExpression, scaleBigMatrix, tmp.dir)
  
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
    nPerm <- sapply(discovery, function(di) {
      modules <- table(moduleAssignments[[di]])
      if (!is.null(excludeModules[[di]])) {
        modules <- modules %sub_nin% excludeModules[[di]]
      }
      if (!is.null(includeModules[[di]])) {
        modules <- modules %sub_in% includeModules[[di]]
      }
      
      nModules <- length(modules)
      # Bonferroni correct for the number of modules of interest in the
      # discovery dataset, multiplied by the number of test datasets, but
      # force a minimum requirement to ensure a reasonable degree of accuracy.
      max(200, requiredPerms(0.05/(nModules*length(test))))
    })
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
  if (.Platform$OS.type == "windows" & nCores > 1) {
    # Quietly load parallel backend packages. Throw our own warning and 
    # continue
    if(suppressWarnings(suppressMessages(require(doParallel)))) {
      # we need an additional thread to monitor and report progress
      if (verbose)  
        nCores <- nCores + 1
      cl <- makeCluster(nCores)
      registerDoParallel(cl)
      on.exit({
        stopCluster(cl)
      }, add=TRUE)
      vCat(verbose, 0, "Running on", nCores - 1, "cores.")
      if ((nCores - 1) > detectCores()) {
        stop(
          "Requested number of threads (", nCores - 1, ") is higher than the ",
          "number of available cores (", detectCores(), "). Using too many ",
          "threads may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end 
      # once the analysis has finished.
      vCat(
        TRUE, 0, file=stderr(),
        "Warning: unable to find 'doParallel' package, running on 1 core.", 
      )
      warning("U")
    }
  } else if (.Platform$OS.type == "unix" & nCores > 1) {
    # Quietly load parallel backend packages. Throw our own warning and 
    # continue
    if(suppressWarnings(suppressMessages(require(doMC)))) {
      # we need an additional thread to monitor and report progress
      if (verbose) 
        nCores <- nCores + 1
      registerDoMC(nCores)
      vCat(verbose, 0, "Running on", nCores - 1, "cores.")
      if ((nCores - 1) > detectCores()) {
        stop(
          "Requested number of threads (", nCores - 1, ") is higher than the ",
          "number of available cores (", detectCores(), "). Using too many ",
          "threads may cause the machine to thrash/freeze."
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
  
  # Since we expect the user to explicitly handle the number of parallel threads,
  # we will disable the potential implicit parallelism on systems where R has
  # been compiled against a multithreaded BLAS, e.g. OpenBLAS. 
  omp_set_num_threads(1)
  blas_set_num_threads(1)
  
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
            di, ", in dataset ", ti, "."
          )
          #---------------------------------------------------------------------
          # Set up variables for this comparison
          #---------------------------------------------------------------------
          
          # Force modules to be character vectors to avoid improper indexing
          if (is.numeric(moduleAssignments[[di]])) { 
            moduleAssignments[[di]] <- as.character(moduleAssignments[[di]])
            names(moduleAssignments[[di]]) <- colnames(coexpression[[di]])
          }
          
          # To simplify later function calls, we need to get a vector of module
          # assignments only for (a) modules of interest and (b) the genes 
          # present in both datasets for those modules.
          overlapGenes <- intersect(
            colnames(coexpression[[di]]), 
            colnames(coexpression[[ti]])
          )
          overlapAssignments <- moduleAssignments[[di]][overlapGenes]
          # Restrict to modules of interest
          modules <- unique(moduleAssignments[[di]])
          if (!is.null(excludeModules[[di]]))
            modules <- modules %sub_nin% excludeModules[[di]]
          if (!is.null(includeModules[[di]]))
            modules <- modules %sub_in% includeModules[[di]]
          overlapAssignments <- overlapAssignments %sub_in% modules 
          overlapModules <- unique(overlapAssignments)
          overlapModules <- overlapModules[orderAsNumeric(overlapModules)]
          
          if (length(overlapAssignments) == 0) {
            warning(
              "No genes for the modules of interest are present in the test",
              " dataset"
            )
            next
          }
          
          # How many genes are present in the test dataset for the modules of 
          # interest?
          genesPres <- table(overlapAssignments)
          modulesWithNoOverlap <- modules %sub_nin% overlapAssignments
          genesPres <- c(genesPres, rep(0, length(modulesWithNoOverlap)))
          names(genesPres)[names(genesPres) == ""] <- modulesWithNoOverlap
          genesPres <- genesPres[orderAsNumeric(names(genesPres))]
          
          # What proportion?
          moduleSizes <- table(moduleAssignments[[di]])
          moduleSizes <- moduleSizes[names(moduleSizes) %sub_in% modules]
          propGenesPres <- genesPres / moduleSizes
          propGenesPres <- propGenesPres[orderAsNumeric(names(propGenesPres))]
          
          # Calculate some basic cross-tabulation statistics so we can assess 
          # which modules in both datasets map to each other, if module
          # detection has also been performed for the test network
          contingency <- NULL
          if (!is.null(moduleAssignments[[ti]])) {
            # Get total number of nodes from each discovery subset in each test subset 
            contingency <- table(
              moduleAssignments[[di]][overlapGenes], 
              moduleAssignments[[ti]][overlapGenes]
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
            testSizes <- table(moduleAssignments[[ti]][overlapGenes])
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
            modGenes <- names(overlapAssignments %sub_in% mi)
            modInds <- match(modGenes, colnames(coexpression[[di]]))
            discProps[[mi]] <- moduleProps(adjacency[[di]], modInds, sge[[di]])
          }
          
          #---------------------------------------------------------------------
          # Calculate the observed value for each test statistic
          #---------------------------------------------------------------------
          observed <- matrix(NA, nrow=nModules, ncol=nStatistics)
          rownames(observed) <- overlapModules
          for (mi in overlapModules) {
            modGenes <- names(overlapAssignments %sub_in% mi)
            discInds <- match(modGenes, colnames(coexpression[[di]]))
            testInds <- match(modGenes, colnames(coexpression[[ti]]))
            testProps <-  moduleProps(adjacency[[ti]], testInds, sge[[ti]])
            stats <- calcStats(
              discProps[[mi]], testProps, 
              coexpression[[di]], discInds,
              coexpression[[ti]], testInds
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
              run.dir <- paste0(".run-progress", getUUID())
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
          foreach(chunk=ichunkTasks(verbose, nPerm[di], nCores)) %maybe_do_par% {
            if (verbose && length(chunk) == 1 && chunk == -1) {
              monitorProgress(nCores - 1, 2, run.dir)
              NULL
            } else {
              if (verbose) {
                conns <- setupParProgressLogs(chunk, nCores - 1, 2, run.dir)
                progressBar <- conns[[1]]
                on.exit(lapply(conns, close))
              } 
              
              # Change attached status of matrices depending on memory requirements
              if (!lowmem) {
                if (!is.null(geneExpression))
                  sge <- lapply(sge, attach.bigMatrix)
                coexpression <- lapply(coexpression, attach.bigMatrix)
                adjacency <- lapply(adjacency, attach.bigMatrix)
                on.exit({
                  if (!is.null(geneExpression))
                    sge <- lapply(sge, detach.bigMatrix)
                  coexpression <- lapply(coexpression, detach.bigMatrix)
                  adjacency <- lapply(adjacency, detach.bigMatrix)
                }, add=TRUE)
              }
              
              #-----------------------------------------------------------------
              # Calculate the module preservation statistics for each module on
              #  a random subset of genes of the same size in the test dataset
              #-----------------------------------------------------------------
              chunkStats <- array(
                NA, dim=c(nModules, nStatistics, length(chunk))
              )
              dimnames(chunkStats)[1:2] <- dimnames(observed)
              dimnames(chunkStats)[[3]] <- paste0("permutation.", chunk)
              for (pi in seq_along(chunk)) {
                for (mi in overlapModules) {
                  modGenes <- names(overlapAssignments %sub_in% mi)
                  discInds <- match(modGenes, colnames(coexpression[[di]]))
                  
                  # Select a random subset of nodes of the same size as the subset 
                  # ss, depending on our null model.
                  modSize <- length(modGenes)
                  if (model == "overlap") {
                    permGenes <- sample(names(overlapGenes), modSize)
                  } else {
                    permGenes <- sample(colnames(coexpression[[ti]]), modSize)
                  }
                  permInds <- match(permGenes, colnames(coexpression[[ti]]))
                  # Ensure crashes aren't fatal
                  tryCatch({
                    permProps <- moduleProps(adjacency[[ti]], permInds, sge[[ti]])
                    chunkStats[mi,,pi] <- calcStats(
                      discProps[[mi]], permProps, 
                      coexpression[[di]], discInds,
                      coexpression[[ti]], permInds
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
              if (colnames(observed)[si] %in% c("mean.adj", "propVarExpl")) {
                order <- FALSE
              } else {
                order <- TRUE
              }
              
              # If the two-sided alternative has been specified, it only makes
              # sense to use it on the correlation-based statistics.
              if (altMatch == 1) {
                if (colnames(observed)[si] %in% c("mean.adj", "propVarExpl")) {
                  alternative <- "greater"
                } else {
                  alternative <- "two.sided"
                }
              } else {
                alternative <- validAlts[altMatch]
              }
              
              # Get the p-values
              p.values[mi, si] <- perm.test(
                nulls[mi, si, ], observed[mi, si], 
                genesPres[mi], length(overlapGenes),
                order=order, alternative=alternative
              )
              
            }
          }
          
          # Order statistics: First density stats, then connectivity
          if (!is.null(sge[[di]])) {
            statOrder <- c(
              "mean.adj", "pve", "cor.coexp", "cor.kIM", "cor.MM",
              "mean.coexp", "mean.MM"
            ) 
          } else {
            statOrder <- c("mean.adj", "cor.kIM", "cor.coexp", "mean.coexp")
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
            nGenesPresent = genesPres,
            propGenesPresent = propGenesPres,
            contingency = contingency
          )
          # remove NULL outputs
          res[[di]][[ti]] <- res[[di]][[ti]][
            !sapply(res[[di]][[ti]], is.null)
          ]
          
          gc()
          on.exit({
            vCat(verbose, 0, "Done!")
          }, add=TRUE)
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
