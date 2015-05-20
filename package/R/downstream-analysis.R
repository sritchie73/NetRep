#' Get the properties of a network module
#' 
#' Calculate the network properties that compose the module preservation 
#' statistics for one or more modules.
#' 
#' @template api_inputs
#' 
#' @return 
#'  A list of network properties for each module of interest.
#'
#' @examples
#' ## Example 1: calculate network properties for a module 
#' ## or gene set.
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
#' # Create bigMatrix objects for each matrix.
#' geA <- as.bigMatrix(geA, "geA_bm")
#' coexpA <- as.bigMatrix(coexpA, "coexpA_bm")
#' adjA <- as.bigMatrix(adjA, "adjA_bm")
#' 
#' # Get the properties of module 2
#' networkProperties(
#'   geA, coexpA, adjA, moduleAssignments, modules="2"
#' )
#' 
#' ## Example 2: calculate the network properties on a user defined 
#' ## subset of the first 10 genes.
#' networkProperties(
#'  geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10]
#' )
#' 
#' ## Example 3: calculate the network properties of an adipose tissue
#' ## module in the liver tissue of the same samples
#' 
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
#' # Store each input type as a list, where each element corresponds
#' # to one of the datasets
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
#' # Get the network properties in the liver tissue for modules 
#' # 3 and 7, which were discovered in the adipose tissue. 
#' networkProperties(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' 
#' @rdname networkProperties
#' @export
networkProperties <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1
) {
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structure
  if (missing(moduleAssignments))
    modules <- "1"
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # Sanity check input for consistency.
  checkSets(
    geneExpression, coexpression, adjacency, moduleAssignments, discovery, test
  )
  
  if (is.null(moduleAssignments[[discovery]]))
    stop("no module assignments in the discovery dataset")
  
  # Temporarily create scaled gene expression set for the calculation of the
  # summary expression profile
  sge <- NULL
  if (!is.null(geneExpression)) {
    tryCatch({
      checkFinite(geneExpression[[test]]) 
    }, error = function(e) {
      stop("Non-finite values encountered for the test dataset gene expression")
    })
    sge <- scaleBigMatrix(geneExpression[[test]], tmp.dir)
  }
  
  # Get the properties for each module of interest 
  res <- lapply(modules, function(mod) {
    # Get the row/column indices of the module in the dataset of interest 
    sub <- moduleAssignments[[discovery]][moduleAssignments[[discovery]] == mod]
    modInds <- match(names(sub), rownames(coexpression[[test]]))
    na.inds <- which(is.na(modInds))
    modInds <- na.omit(modInds)
    
    if (length(modInds) == 0) {
      stop(
        "none of the genes for module ", mod, 
        " are present in the test dataset"
      )
    }
    
    tryCatch({
      checkFinite(adjacency[[test]]) 
    }, error = function(e) {
      stop("Non-finite values encountered for the test dataset adjacency")
    })
    
    # Get the properties calculated from the gene expression
    geProps <- NULL
    if (!is.null(sge)) {
      geProps <- dataProps(sge, modInds)
      # rename for clarity
      names(geProps) <- c("summaryExpression", "moduleMembership", "propVarExpl")
      geProps[[2]] <- insert.nas(geProps[[2]], na.inds)
      names(geProps[[1]]) <- rownames(sge)
      names(geProps[[2]]) <- names(sub)
    }
    
    # Get the properties calculated from the adjacency.
    adjProps <- adjProps(adjacency[[test]], modInds)
    names(adjProps) <- c("connectivity", "density")
    adjProps[[1]] <- insert.nas(adjProps[[1]], na.inds)
    names(adjProps[[1]]) <- names(sub)
    
    c(geProps, adjProps)
  })
  if (length(res) == 1) {
    res <- res[[1]]
  } else {
    names(res) <- modules
  } 
  res
}

#' Order genes and modules within a network.
#' 
#' Get the ordering of genes by intramodular connectivity and ordering of 
#' modules by the similarity of their summary expression profiles. This function
#' is used in the \code{\link[=plotNetwork]{plotting functions}} but may also be
#' useful for other types of downstream analysis: gene connectivity within a
#' module is strongly correlated with a gene's biological importance to that
#' module \emph{(1)}.
#' 
#' @template api_inputs
#' 
#' @references
#' \enumerate{
#'    \item{
#'      Langfelder, P., Mischel, P. S. & Horvath, S. \emph{When is hub gene 
#'      selection better than standard meta-analysis?} PLoS One \strong{8}, 
#'      e61505 (2013).
#'    }
#' }
#' 
#' @return
#'  A vector of ordered gene names (descending order).
#'  
#' @examples
#' ## Example 1: get the ordering of genes for a single module
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
#' # Create bigMatrix objects for each matrix.
#' geA <- as.bigMatrix(geA, "geA_bm")
#' coexpA <- as.bigMatrix(coexpA, "coexpA_bm")
#' adjA <- as.bigMatrix(adjA, "adjA_bm")
#' 
#' # Get the order of genes in module 2
#' orderGenes(
#'   geA, coexpA, adjA, moduleAssignments, modules="2"
#' )
#' 
#' ## Example 2: get the order of genes of an arbitrary subset
#' ## (the first 10 genes)
#' orderGenes(
#'  geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10]
#' )
#' 
#' ## Example 3: get the ordering of genes for two adipose 
#' ## tissue modules in the liver tissue of the same samples
#' 
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
#' # Store each input type as a list, where each element corresponds
#' # to one of the datasets
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
#' # Get the order of genes in the liver tissue for modules 
#' # 3 and 7, which were discovered in the adipose tissue. 
#' orderGenes(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' 
#' @rdname orderGenes
#' @export
orderGenes <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1
) {
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules,
    discovery, test
  )
  if (missing(modules) || length(modules) == 1) {
    props <- list(props)
  }
  
  # order modules
  moduleOrder <- 1
  if (length(props) > 1) {
    if (!is.null(geneExpression[[test]])) {
      # Create a matrix of summary expression profiles to measure the similarity
      seps <- matrix(
        0, ncol=length(props), nrow=length(props[[1]]$summaryExpression)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summaryExpression
      }
      moduleOrder <- hclust(as.dist(1-cor(seps)))$order
    } else {
      warning(
        "No gene expression provided, modules will be ordered as provided"
      )
      moduleOrder <- seq_along(props)
    }
  } 
  
  # order genes
  foreach(mi = moduleOrder, .combine=c) %do% {
    names(sort(props[[mi]]$connectivity, decreasing=TRUE))
  }
}