#' Assessing Module Preservation
#' 
#' @description
#' This is the main package function that performs the module preservation test.
#' It requires well formed input: please see details for instructions on input.
#' 
#' In future releases, more user friendly interfaces will be provided and this
#' function will not be exported to the namespace
#' 
#' @details
#' \subsection{Input Instructions:}{
#'  Any function argument ending in \code{"Sets"} expects a list, where each 
#'  element corresponds to a dataset. The order of datasets should match across
#'  all \code{"Sets"} arguments. 
#'  
#'  \code{exprSets}, \code{coexpSets}, and \code{adjSets} should contain the 
#'  file paths to the file-backed \code{\link[bigmemory]{big.matrix}} objects
#'  for the gene expression, gene coexpression, and gene adjacencies matrices,
#'  respectively, for each dataset. Matrices in R can be converted to 
#'  \code{big.matrix} objects via \code{\link[bigmemory]{as.big.matrix}}, or 
#'  alternatively read in from disk as via \code{\link[bigmemory]{read.big.matrix}}.
#'  After the matrices have been converted to \code{big.matrix} objects the 
#'  original matrices or resulting pointers returned by \code{as.big.matrix} and
#'  \code{read.big.matrix} should be removed via \code{\link{rm}}, and the 
#'  garbage collector invoked \code{\link{gc}} to prevent unnecessary memory 
#'  usage from the left over objects being passed around.
#'  
#'  A nice benefit of creating file-backed \code{big.matrix} objects is that 
#'  they can be attached to any future R session with zero loading time.
#'  
#'  When converting to \code{big.matrix} objects, it is important to ensure that
#'  the column and rownames are ignored to optimise speed. When converting 
#'  matrices in the interactive R session via \code{as.big.matrix}, these should
#'  be first removed by setting the \code{\link{dimnames}} to \code{NULL}. When
#'  reading in the objects from a file using \code{read.big.matrix} the colnames
#'  can be ignored by specifying \code{skip = 1}, and rownames can be ignored 
#'  through a combination of \code{has.row.name = TRUE} and 
#'  \code{ignore.row.names = TRUE}.
#'  
#'  The  \code{geneNameSets} argument is instead used
#'  to hold this information. \code{netRepMain} therefore expects that for any
#'  given dataset the order of genes is the same across the expression data,
#'  coexpression, and gene adjacency matrices. The gene expression matrices are
#'  expected to be in sample x gene format; that is, the genes are the columns.
#'  Each element of \code{geneNameSets} should contain a vector of unique 
#'  identifiers for each probe/gene in the corresponding dataset. Although 
#'  \code{netRepMain} does not require any information about the sample order,
#'  it may be useful to store the rownames of the gene expression data in a 
#'  similar variable, e.g. \code{sampleNameSets}, for downstream analysis.
#'  
#'  For each dataset where module discovery has been performed (hereby referred 
#'  to as the \emph{discovery} dataset), \code{geneModuleSets} should contain a
#'  vector assigning each probe to a module (i.e. the names of the vector are 
#'  the same as \code{geneNameSets}). The entries of \code{geneModuleSets}
#'  corresponding to the datasets the module preservation is being tested in
#'  (hereby referred to as the \emph{test} dataset) should be \code{NULL}. If
#'  module detection has been performed in the \emph{test} dataset as well, then
#'  providing the module assignments will allow \code{netRepMain} to return
#'  a contingency table of module overlap between the two datasets.
#'  
#'  For datasets where only a specific module, or subset of modules, are of
#'  interest, \code{includeSets} should contain a vector of those module labels.
#'  Similarly, \code{excludeSets} can be used to ignore one or more modules.
#' }
#' 
#' \subsection{Parallelism:}{
#'  If \code{verbose} is \code{TRUE} and a parallel backend has been registered,
#'  then one of the cores will be reserved for reporting the function's progress.
#'  In this case an additional core should be registered, even if it means 
#'  registering more "cores" than physically exist: the core reporting the 
#'  progress will have 0 memory and CPU overhead. Note that if only 2 cores have
#'  been registered, this is the same as running on a single core.
#'  
#'  \code{\link{combineNulls}} may be useful for combining the resulting null
#'  distributions from multiple runs (e.g. when parallelising on multiple 
#'  machines without shared memory).
#' }
#' 
#' \subsection{Memory Requirements:}{
#'  \code{netRep} operates on \code{\link[bigmemory]{big.matrix}} objects. This
#'  has two nice consequences:
#'  \enumerate{
#'    \item{
#'      R does not need to load all matrices into RAM, and the individual
#'      matrices in \code{exprSets}, \code{coexpSets}, and \code{adjSets} may also be 
#'      much larger than RAM.
#'    }
#'    \item{
#'      The matrices are accessible from shared memory, meaning that 
#'      \code{netRep} may be run in parallel without having to copy matrices to
#'      each spawned R session, and without storing multiple copies of the 
#'      matrices being operated on.
#'    }
#'  }
#'  
#'  Total runtime is dictated by the available RAM. Due to the nature of 
#'  permutation procedures, the overall runtime for a single core will be 
#'  minimised if all relevant matrices can be stored in RAM for each core.
#'  
#'  When the \code{lowmem} argument is set to \code{TRUE}, \code{netRep} will 
#'  instead consume the minimum possible RAM required for each permutation. This
#'  makes calculations roughly two and a half times slower, but consumes at least
#'  three times less RAM (at peak RAM usage). It is therefore recommended to set
#'  \code{lowmem} to \code{TRUE} on machines with limited RAM, as it will make
#'  parallelisation over many more cores possible without thrashing the machine.
#' }
#' 
#' \subsection{Module Preservation Statistics:}{
#' There are seven module preservation statistics \emph{(1)}:
#'  \enumerate{
#'    \item{\code{mean.adj}:}{
#'      The mean adjacency, or module density, measures how densely connected a
#'      module is in the test dataset.
#'    }
#'    \item{\code{pve}:}{
#'      The pve, short for "the proportion of variance explained in the 
#'      underlying gene expression data for the module by its summary
#'      expression profile in the \emph{test} dataset". The summary expression 
#'      profile is calculated as the first eigenvector from a principal
#'      component analysis on the module's (scaled) gene expression data.
#'      
#'      The summary expression profile is commonly referred to as the "module
#'      eigengene (ME)", and \code{pve} has previously been abbreviated as 
#'      \code{propVarExpl} \emph{(1)}.
#'    }
#'    \item{\code{cor.coexp}:}{
#'      The correlation of coexpression measures how correlated the coexpression 
#'      patterns are across both datasets. It has previously been referred to as
#'      the "correlation of correlation" and abbreviated as \code{cor.cor} 
#'      \emph{(1)}.
#'    }
#'    \item{\code{cor.kIM}:}{
#'      The correlation of intramodular connectivity, measures how correlated
#'      the intramodular connectivity is across datasets. The intramodular 
#'      connectivity is quantified as the sum of adjacencies for a gene to all
#'      other genes in the module.
#'    }
#'    \item{\code{cor.MM}:}{
#'      The correlation of intramodular module membership measures how 
#'      correlated the module membership is across datasets. The module 
#'      membership is quantified as the correlation between each gene and the
#'      summary expression profile. The module membership and summary expression
#'      profile are calculated independently in each dataset. It has previously
#'      been abbreviated as \code{cor.kME} \emph{(1)}.
#'    }
#'    \item{\code{mean.coexp}:}{
#'      The mean sign-aware coexpression measures the average 
#'      coexpression for a module in the \emph{test} dataset, multiplied by the
#'      sign of the coexpression in the \emph{discovery dataset}. It has 
#'      previously been referred to as the "mean sign-aware correlation" and
#'      abbreviated as \code{mean.cor} \emph{(1)}.
#'    }
#'    \item{\code{mean.MM}:}{
#'      The mean sign-aware module membership measures the average module 
#'      membership in the \emph{test} dataset, multiplied by the sign of the 
#'      module membership in the \emph{discovery} dataset. It has previously
#'      been abbreivated as \code{mean.kME}.
#'    }
#'  }
#' }
#' 
#' \subsection{Statistical Tests:}{
#'  When a module is preserved, one would expect the module preservation 
#'  statistics to be larger than expected by chance \emph{(1)}. This is the 
#'  default behaviour. 
#'  
#'  \code{\link{requiredPower}} has been provided to determine the minimum 
#'  number of permutations required to assess module preservation at any given
#'  signficance threshold. Due to the nature of permutation testing, this is 
#'  also the \emph{resolution} of the permutation test; both the smallest 
#'  p-value obtainable, and the increment size between the range of possible
#'  p-values. In general it is advisable to run more permutations than 
#'  necessary for any given significance threshold. 
#'  
#'  There are two interpretations of the null hypothesis: first, only the probes
#'  or genes present in both datasets are considered. In this case probes 
#'  present in the test dataset but not in the discovery are considered "missing
#'  information", so are excluded from the random sampling (\code{null =
#'  "overlap"}). The second is that the module preservation statistics should be
#'  more extreme than any random subset of the test dataset, even when
#'  considering probes not present in the discovery dataset (\code{null =
#'  "all"}). The former is the default behaviour, but the latter may be useful 
#'  when performing multiple species comparisons. 
#'  
#'  One may also be interested in examining cases where the module preservation
#'  statistics are much smaller than expected by chance. These also imply that 
#'  the module is structurally distinct from the rest of the
#'  dataset. Extremely small values for the \emph{pve} or \emph{mean.adj}
#'  indicate that there is less concordance amongst the genes than expected by
#'  chance. Extreme negative values for the \emph{cor.kIM}, \emph{cor.MM}, 
#'  \emph{cor.coexp}, \emph{mean.MM}, or \emph{mean.coexp} indicate that the 
#'  gene-gene relationships are reversed across datasets. Extreme negative 
#'  values for the \emph{cor.coexp} and \emph{mean.coexp} should occur in
#'  cases where the coexpression has flipped: i.e. positively correlated genes
#'  are now negatively correlated and vice versa. Extreme negative values for 
#'  \emph{cor.kIM} indicate that the order of gene importance (ranked by 
#'  intramodular connectivity) is reverse. Extreme negative values for
#'  the \emph{mean.MM} and \emph{cor.MM} may indicate that the summary 
#'  expression is being dominated by a few genes that are negatively
#'  correlated with the rest of the module in one of the datasets (although
#'  one should excercise caution here, as the summary expression profile may
#'  not be representitive of the gene expression if the \emph{pve} is low).
#'  
#'  If one wishes to test for cases where the module preservation statistics are
#'  smaller than expected by chance, \code{alternative} should be set to "less".
#'  A two-sided test may also be performed (\code{alternative = "two.sided"}),
#'  but this requires twice as many permutations to detect extreme observations
#'  at any given threshold.
#'  
#'  P-values may be calculated post-hoc using \code{\link{perm.test}} if one
#'  wishes to evaluate an alternative hypothesis later on the returned null
#'  distributions.
#'}
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
#'   }
#' 
#' @param exprSets \code{NULL}, or a list of \code{\link{big.matrix}}: one for 
#'   each dataset (See Input Instructions in Details).
#' @param coexpSets a list of \code{\link{big.matrix}} objects: one for each 
#'   dataset, corresponding to the pairwise-gene correlation (See Input
#'   Instructions in Details).
#' @param adjSets a list of \code{\link{big.matrix}} objects: one for each
#'   dataset, corresponding to the adjacency matrix of edge weights between each
#'   pair of nodes in the network (See Input Instructions in Details).
#' @param geneNameSets \code{NULL}, or a list of character vectors giving the 
#'   node names for the corresponding element of \code{exprSets}, \code{coexpSets},
#'   and \code{adjSets} (See Input Instructions in Details).
#' @param geneModuleSets a list, whose elements are \code{NULL} for each 
#'   \emph{test} dataset, and a vector for each \emph{discovery} dataset 
#'   assigning each node to a sub-network/cluster/module/component (See Input
#'   Instructions in Details).
#' @param discovery a numeric vector indicating which elements of
#'   \code{exprSets} and/or \code{adjSets} are to be treated as the 
#'   \emph{discovery} datasets.
#' @param test a numeric vector indicating which elements of \code{exprSets}
#'   and/or \code{adjSets} are to be treated as the \emph{test} datasets.
#' @param nPerm either a numeric vector or a list of numeric vectors specifying
#'   the number of permutations to use for all datasets, or each discovery 
#'   dataset. If a list of numeric vectors is supplied it should be in the same
#'   format as the other "Sets" arguments, see Input Instructions in Details.
#' @param excludeSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   skip.
#' @param includeSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   include.
#' @param null the type of null model, either "overlap" or "all". See section
#'  in Details on Statistical Tests.
#' @param alternative The type of module preservation test to perform. Must be 
#'   one of "greater" (default), "less" or "two.sided". See section in details
#'   on Statistical Tests.
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#' @param indent numeric; a positive value indicating the indent level to start
#'   the output at. Defaults to 0. Each indent level adds two spaces to the 
#'   start of each line of output.
#' @param simplify logical; if \code{TRUE}, simplify the structure of the output
#'  list if possible (see Return Value).
#' @param lowmem logical; should memory usage be minimised? Useful on machines
#'  with limited RAM when running in parallel. See section on Memory Usage in
#'  Details.
#'  
#' @return
#'  The returned data structure is organised as a nested list of lists, which 
#'  should be accessed as \code{results[[discovery]][[test]]}. These are 
#'  typically numeric indices, i.e. if the \emph{discovery} dataset is the first in 
#'  the \code{"Sets"} list (see Input Instructions in details) and the 
#'  \emph{test} dataset is the second, then the results will be stored in 
#'  \code{results[[1]][[2]]}, while \code{results[[2]][[1]]} will be empty, and 
#'  so on.
#'  
#'  If \code{simplify} is set to \code{TRUE}, then this structure will be
#'  simplified as much as possible depending on the combination of dataset
#'  comparisons that have been performed. 
#'  
#'  For each dataset-comparison a list of the following objects are returned:
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
#'    \item{\code{Z.scores}:}{
#'      A matrix of Z scores for the \code{observed} module preservation
#'      statistics.
#'    }
#'    \item{\code{genesPresent}:}{
#'      A vector containing the number of genes that are present in the test
#'      dataset for each module.
#'    }
#'    \item{\code{propGenesPresent}:}{
#'      A vector containing the proportion of genes present in the test dataset
#'      for each module. Modules where this is less than 1 should be 
#'      investigated further before making judgements about preservation to 
#'      ensure that the missing genes are not the most connected ones.
#'    }
#'    \item{\code{contigency}:}{
#'      If the corresponding \code{geneModuleSets} for the \emph{test} dataset
#'      has also been supplied, then a contigency table is returned, showing
#'      the overlap between modules in each dataset. Rows correspond to modules
#'      in the \emph{discovery} dataset, columns to modules in the \emph{test}
#'      dataset.
#'    }
#'  }
#'
#' @import foreach
#' @import RhpcBLASctl
#' @export
netRepMain <- function(
  exprSets=NULL, coexpSets, adjSets, geneNameSets=NULL, 
  geneModuleSets, discovery, test, nPerm=10000, excludeSets=NULL, 
  includeSets=NULL, null="overlap", alternative="greater",
  verbose=TRUE, indent=0, simplify=TRUE, lowmem=FALSE
) {
  # The following declarations are for iterators declared inside each foreach 
  # loop. Declarations are required to satisfy NOTES generated by R CMD check, 
  # and also serve as useful documentation for the reader of the source code.
  di <- NULL # discovery dataset index in the "Sets" lists.
  ti <- NULL # test dataset index in the "Sets" lists.
  ii <- NULL # iterator over the subsets
  jj <- NULL # iterator over the statistics
  kk <- NULL # iterator over the permutations
  
  # Identify the null model to use 
  nullModels <- c("overlap", "all")
  if (is.na(pmatch(null, nullModels))) {
    stop("overlap must match one of ", paste(nullModels, collapse=" "))
  }
  model <- pmatch(null, nullModels)
  
  # Determine the type of test we're performing
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  
  # Determine set up of worker nodes (if applicable)
  nCores <- getDoParWorkers()
  if (verbose & nCores > 1) {
    vCat(verbose, indent, "Verbose output is TRUE; reserving one core as the", 
         "master, to report the progress of the worker cores who do the bulk",
         "of the calculation.")
    nWorkers <- nCores - 1
  } else {
    nWorkers <- nCores
  }
  vCat(verbose, indent, "Running with", nWorkers, "worker cores.")
  
  # Since we expect the user to explicitly handle the number of parallel threads,
  # we will disable the potential implicit parallelism on systems where R has
  # been compiled against a multithreaded BLAS, e.g. OpenBLAS. 
  omp_set_num_threads(1)
  blas_set_num_threads(1)
  
  # Note to self: if we do any linear algebra on the full networks in 
  # preparation at a later date, we'll first need to set the number of threads
  # to `nCores`.
  
  nNets <- length(geneModuleSets)
  
  # Set up names for more informative output, if the datasets and their 
  # corresponding networks have been given names.
  if (is.null(names(geneModuleSets))) {
    setNames <- 1:nNets
  } else {
    setNames <- names(geneModuleSets)
  }
  
  # Set up temporary directory for working big.matrix objects, and permutation
  # results.
  while (TRUE) {
    obj.dir <- paste0(".temp-objects", getUUID())
    # Handle the infintesimally small chance of a UUID collision
    tryCatch({
      dir.create(obj.dir)
      break
    }, warning = function(w) {
      if(!grepl("already exists", w$message)) {
        break
      }
    })
  }  
  on.exit({
    unlink(obj.dir, recursive=TRUE)
  }, add=TRUE)
  
  scaledSets <- NULL
  if (!is.null(exprSets)) {
    scaledSets <- rep(list(NULL), length(exprSets))
  }
  
  # Set up return list 
  res <- rep(list(NULL), nNets)
  res <- lapply(res, function(x) { 
    l <- rep(list(NULL), nNets)
    names(l) <- names(geneModuleSets)
    l
  })
  names(res) <- names(geneModuleSets)
  
  # Iterate pairwise over datasets, comparing those marked "discovery"
  # with each marked as "test".
  for (di in seq_len(nNets)) {
    for (ti in seq_len(nNets)) {
      if ((di %in% discovery) & (ti %in% test) & (di != ti)) {
        tryCatch({
          # Set up return list
          res[[di]][[ti]] <- rep(list(NULL), 5)
          
          hasDat <- FALSE # Flag useful for later
          
          # Output messages
          vCat(verbose, indent, sep="", 
               "Calculating preservation of network subsets from dataset ",
               setNames[di], ", in dataset ", setNames[ti], ".")
          
          # Check input is ok
          if (length(nPerm) > 1) {
            if(!(class(nPerm[[di]]) %in% c("numeric", "integer", "double"))) {
              stop("'nPerm' for dataset", di, "has not been supplied")
            }
            thisPerm <- nPerm[[di]]
          } else {
            thisPerm <- nPerm
          }
          
          vCat(verbose, indent+1, "Checking matrices...")
          checkFinite(coexpSets[[di]])
          checkFinite(adjSets[[di]])
          checkFinite(coexpSets[[ti]])
          checkFinite(adjSets[[ti]])
          if (!is.null(exprSets[[di]]) && !is.null(exprSets[[ti]])) {
            checkFinite(exprSets[[di]])
            checkFinite(exprSets[[ti]])
            hasDat <- TRUE
          } 
          
          # Create scaled data 
          if (hasDat) {
            vCat(verbose, indent+1, "Creating temporary scaled gene expression...")
            if (is.null(scaledSets[[di]])) {
              scaledSets[[di]] <- scaleBigMatrix(exprSets[[di]], obj.dir)
            } 
            if (is.null(scaledSets[[ti]])) {
              scaledSets[[ti]] <- scaleBigMatrix(exprSets[[ti]], obj.dir)
            }
          }
          
          if (lowmem) {
            scaledDisc <- scaledSets[[di]]
            scaledTest <- scaledSets[[ti]]
            corDisc <- coexpSets[[di]]
            corTest <- coexpSets[[ti]]
            adjDisc <- adjSets[[di]]
            adjTest <- adjSets[[ti]]
          } else {
            if (hasDat) {
              scaledDisc <- attach.big.matrix(scaledSets[[di]])
              scaledTest <- attach.big.matrix(scaledSets[[ti]])
            } else {
              scaledDisc <- NULL
              scaledTest <- NULL
            }
            corDisc <- attach.big.matrix(coexpSets[[di]])
            corTest <- attach.big.matrix(coexpSets[[ti]])
            adjDisc <- attach.big.matrix(adjSets[[di]])
            adjTest <- attach.big.matrix(adjSets[[ti]])
          }
          
          # Get a vector of nodes which are present in both datasets. Depends on 
          # the combination of data input provided.
          vCat(verbose, indent+1, "Extracting information about node overlap...")
          oNodes <- geneNameSets[[di]] %sub_in% geneNameSets[[ti]]
          tNodes <- geneNameSets[[ti]]
          
          # Force the node indices to be strings, even if the provided identifiers
          # are integers.
          oNodes <- as.character(oNodes)
          tNodes <- as.character(tNodes)
          
          if (length(oNodes) == 0) {
            warning("No nodes in dataset ", setNames[di],  
                    " are present in dataset ", setNames[ti], ", skipping.")
            next
          }
          
          # Compute information about the network subsets, their size, and what 
          # proportion of each is overlapping the test network.
          dSubsets <- unique(geneModuleSets[[di]])
          oSubsets <- unique(geneModuleSets[[di]][oNodes])
          # Only look at subsets of interest, if specified:
          if (!is.null(excludeSets[[di]])) {
            dSubsets <- dSubsets %sub_nin% excludeSets[[di]]
            oSubsets <- oSubsets %sub_nin% excludeSets[[di]]
          }
          if (!is.null(includeSets[[di]])) {
            dSubsets <- dSubsets %sub_in% includeSets[[di]]
            oSubsets <- oSubsets %sub_in% includeSets[[di]]
          }
          # Get the size of each of the subsets of interest in the discovery dataset
          dSizes <- table(geneModuleSets[[di]])
          dSizes <- dSizes[names(dSizes) %in% dSubsets] 
          # Get the size of the overlap of each of the subsets of interest between
          # the discovery and test datasets.
          genesPres <- table(geneModuleSets[[di]][oNodes])
          genesPres <- genesPres[names(genesPres) %in% oSubsets]
          # Need to handle the case where a subset of interest has no overlap.
          genesPres <- c(genesPres, rep(0, length(dSubsets %sub_nin% oSubsets)))
          names(genesPres) <- c(names(genesPres) %sub_nin% "", dSubsets %sub_nin% oSubsets)
          propGenesPres <- genesPres[names(dSizes)]/dSizes
          
          # How many network properites and statistics will we return? 
          # Numbers required for data structure allocation
          nStats <- ifelse(hasDat, 7, 4)
          nSubsets <- length(oSubsets)
          
          # Calculate some basic cross-tabulation statistics so we can assess
          # which modules in both datasets map to each other, if module detection
          # has also been performed for the test network
          if (!is.null(geneModuleSets[[ti]])) {
            # Get total number of nodes from each discovery subset in each test subset 
            subsetOverlap <- table(
              geneModuleSets[[di]][oNodes], 
              geneModuleSets[[ti]][oNodes]
            )
            # filter on subsets the user cares about
            subsetOverlap <- subsetOverlap[dSubsets,]
            # order the tables
            tryCatch({
              # For modules that are integer coded, make sure they're numerically
              # ordered, not alphabetically.
              rOrder <- order(as.integer(rownames(subsetOverlap)))
              cOrder <- order(as.integer(colnames(subsetOverlap)))
            }, warning = function(w) {
              # If we can't cast to an integer, sort normally.
              rOrder <- order(rownames(subsetOverlap))
              cOrder <- order(colnames(subsetOverlap))
            })
            subsetOverlap <- subsetOverlap[rOrder, cOrder]
            # add information about sizes of both the discovery and test subsets
            tryCatch({
              # For modules that are integer coded, make sure they're numerically
              # ordered, not alphabetically.
              dSubSizes <- table(geneModuleSets[[di]][oNodes])
              dSubSizes <- dSubSizes[order(as.integer(names(dSubSizes)))]
            }, warning=function(w) { 
              # ignore
            })
            tryCatch({
              # For modules that are integer coded, make sure they're numerically
              # ordered, not alphabetically.
              tSubSizes <- table(geneModuleSets[[ti]][oNodes])
              tSubSizes <- tSubSizes[order(as.integer(names(tSubSizes)))]
            }, warning=function(w) { 
              # ignore
            })

            subsetOverlap <- cbind(dSubSizes, subsetOverlap)
            subsetOverlap <- rbind(c(NA, tSubSizes), subsetOverlap)
            rownames(subsetOverlap)[1] <- "size"
            colnames(subsetOverlap)[1] <- "size"
          } else {
            subsetOverlap <- NULL
          }
          
          # Obtain the topological properties for each network subset in the
          # discovery dataset, we only want to calculate these once!
          vCat(verbose, indent+1, "Calculating observed test statistics...")
          discProps <- rep(list(NULL), nSubsets)
          for (ii in seq_along(oSubsets)) {
            sNodes <- names(which(geneModuleSets[[di]][oNodes] == oSubsets[ii]))
            # get the indices in the underlying data and adjacency matrices for 
            # the subset nodes.
            subsetInd <- match(sNodes, geneNameSets[[di]])
            
            props <- moduleProps(adjDisc, subsetInd, scaledDisc, lowmem)
            discProps[[ii]] <- props
          }
          names(discProps) <- oSubsets
          
          # Now calculate the observed value for each network statistic
          observed <- matrix(NA, nrow=nSubsets, ncol=nStats)
          rownames(observed) <- oSubsets
          for (ii in seq_along(oSubsets)) {
            ss <- as.character(oSubsets[ii])
            sNodes <- names(which(geneModuleSets[[di]][oNodes] == ss))
            testInd <- match(sNodes, geneNameSets[[ti]])
            discInd <- match(sNodes, geneNameSets[[di]])
            testProps <- moduleProps(adjTest, testInd, scaledTest, lowmem)
            stats <- c(
              calcSplitTestStats(discProps[[ss]], testProps),
              calcSharedTestStats(corDisc, discInd, corTest, testInd, lowmem)
            )
            observed[ii,] <- stats
            colnames(observed) <- names(stats)
          }
          
          # Clean up master process and free memory
          rm(scaledDisc, scaledTest, corDisc, corTest, adjDisc, adjTest)
          gc()
          
          # Calculate the null distribution for each of the statistics.
          vCat(
            verbose, indent+1, "Calculating null distributions with", thisPerm, 
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
          foreach(chunk=ichunkTasks(verbose, thisPerm, nCores)) %maybe_do_par% {
            if (verbose & length(chunk) == 1) {
              if (chunk == -1) {
                monitorProgress(nWorkers, indent+2, run.dir)
                NULL
              }
            } else {
              if (verbose) {
                conns <- setupParProgressLogs(chunk, nWorkers, indent+2, run.dir)
                progressBar <- conns[[1]]
                on.exit(lapply(conns, close))
              } 
              
              if (lowmem) {
                scaledTest <- scaledSets[[ti]]
                corDisc <- coexpSets[[di]]
                corTest <- coexpSets[[ti]]
                adjTest <- adjSets[[ti]]
              } else {
                if (hasDat) {
                  scaledTest <- attach.big.matrix(scaledSets[[ti]])
                } else {
                  scaledTest <- NULL
                }
                corDisc <- attach.big.matrix(coexpSets[[di]])
                corTest <- attach.big.matrix(coexpSets[[ti]])
                adjTest <- attach.big.matrix(adjSets[[ti]])
              }
              
              chunkStats <- array(NA, dim=c(nSubsets, nStats, length(chunk)))
              dimnames(chunkStats)[1:2] <- dimnames(observed)
              dimnames(chunkStats)[[3]] <- paste0("permutation.", chunk)
              for (kk in seq_along(chunk)) {
                for (ii in seq_along(oSubsets)) {
                  ss <- as.character(oSubsets[ii])
                  # Select a random subset of nodes of the same size as the subset 
                  # ss, depending on our null model.
                  if (model == "overlap") {
                    permNames <- sample(oNodes, size=genesPres[ss])
                  } else {
                    permNames <- sample(tNodes, size=genesPres[ss])
                  }
                  
                  permInd <- match(permNames, geneNameSets[[ti]])
                  sNodes <- names(which(geneModuleSets[[di]][oNodes] == ss))
                  discInd <- match(sNodes, geneNameSets[[di]])
                  
                  tryCatch({ 
                    testProps <- moduleProps(
                      adjTest, permInd, scaledTest, lowmem
                    )
                    chunkStats[ii,,kk] <- calcStats(
                      discProps[[ss]], testProps, corDisc, discInd, corTest, 
                      permInd, lowmem
                    )
                    rm(testProps)
                    gc()
                  }, error = function(e) {
                    warning(
                      "Calculation for subset ", oSubsets[ii], " failed on ",
                      "permutation ", chunk[kk], " with error message:\n",
                      e$message
                    )
                  })
                }
                # Update the progress at the end of the loop.
                if (verbose) {
                  updateParProgress(progressBar, chunk[kk])
                  if (nCores == 1) {
                    reportProgress(indent+2, run.dir)
                    if (chunk[kk] == thisPerm) {
                      cat("\n")
                    }
                  }
                }
              }
              chunkNum <- ceiling(chunk[1]/length(chunk))
              permFile <- paste0("chunk", chunkNum, "permutations.rds")
              saveRDS(chunkStats, file.path(obj.dir, permFile))
            }
          }
          # Get rid of discProps, we no longer need this!
          rm(discProps)
          gc()
          
          # Load in results
          nulls <- array(NA, dim=c(nSubsets, nStats, thisPerm))
          dimnames(nulls)[[3]] <- rep("", dim(nulls)[3])
          chunkFiles <- list.files(obj.dir, "chunk[0-9]*permutations.rds")
          offset <- 1
          for (cf in chunkFiles) {
            chunk <- readRDS(file.path(obj.dir, cf))
            nCPerm <- dim(chunk)[3]
            nulls[,,offset:(offset+nCPerm-1)] <- chunk
            dimnames(nulls)[1:2] <- dimnames(chunk)[1:2]
            dimnames(nulls)[[3]][offset:(offset+nCPerm-1)] <- dimnames(chunk)[[3]]
            offset <- offset + nCPerm
          }        
          
          # Calculate the p-value for the observed statistic based on the null 
          # distribution
          vCat(verbose, indent+1, "Calculating P-values...")
          p.values <- matrix(NA, nrow=nSubsets, ncol=nStats)
          dimnames(p.values) <- dimnames(observed)
          Z.scores <- p.values
          for (ii in seq_along(oSubsets)) {
            for (jj in seq_len(nStats)) {
              
              # Does the order of nodes in each permutation affect the statistic?
              if (colnames(observed)[jj] %in% c("mean.adj", "propVarExpl")) {
                order <- FALSE
              } else {
                order <- TRUE
              }
              
              # If the two-sided alternative has been specified, it only makes
              # sense to use it on the correlation-based statistics.
              if (altMatch == 1) {
                if (colnames(observed)[jj] %in% c("mean.adj", "propVarExpl")) {
                  alternative <- "greater"
                } else {
                  alternative <- "two.sided"
                }
              } else {
                alternative <- validAlts[altMatch]
              }
              
              # Get the p-values
              p.values[ii, jj] <- perm.test(
                nulls[ii, jj, ], observed[ii, jj], 
                genesPres[rownames(p.values)[ii]], length(oNodes),
                order=order, alternative=alternative
              )
              
              Z.scores[ii, jj] <- observed[ii, jj]
            }
          }
          Z.scores <- Z.scores - apply(nulls, c(1,2), mean)
          Z.scores <- Z.scores/apply(nulls, c(1,2), sd)
                    
          # Collate results
          # First order output nicely
          tryCatch({
            # For modules that are integer coded, make sure they're numerically
            # ordered, not alphabetically.
            arrOrder <- order(as.integer(rownames(nulls)))
            vOrder <- order(as.integer(names(propGenesPres)))
          }, warning = function(w) {
            # If we can't cast to an integer, sort normally.
            arrOrder <- order(rownames(nulls))
            vOrder <- order(propGenesPres)
          })

          # Order statistics: First density stats, then connectivity
          if (hasDat) {
            statOrder <- c(
              "mean.adj", "pve", "cor.coexp", "cor.kIM", "cor.MM",
              "mean.coexp", "mean.MM"
            ) 
          } else {
            statOrder <- c("mean.adj", "cor.kIM", "cor.coexp", "mean.coexp")
          }
          
          vCat(verbose, indent+1, "Collating results...")
          
          # Collate results
          res[[di]][[ti]] <- list(
            observed = observed[arrOrder, statOrder],
            nulls = nulls[arrOrder, statOrder,],
            p.values = p.values[arrOrder, statOrder],
            Z.scores = Z.scores[arrOrder, statOrder],
            genesPresent = genesPres[vOrder],
            propGenesPresent = propGenesPres[vOrder]
          )
          
          if(!is.null(subsetOverlap)) {
            res[[di]][[ti]][[length(res[[di]][[ti]]) + 1]] <- subsetOverlap
            names(res[[di]][[ti]])[length(res[[di]][[ti]])] <- "contingency"
          }
                    
          vCat(verbose, indent+1, "Cleaning up temporary objects...")
          unlink(run.dir, recursive=TRUE)
          gc()
          vCat(verbose, indent, "Done!")
        }, error=function(e) {
          warning(
            "Failed with error:\n", e$message, "\nSkipping to next comparison"
          )
        })
      }
    }
  }
  if (simplify) {
    # remove doubly-nested list structure where possible
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

