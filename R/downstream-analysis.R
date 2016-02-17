#' Calculate the topological properties for a network module
#' 
#' Calculates the network properties used to assess module preservation for one
#' or more modules in a user specified dataset.
#' 
#' @inheritParams common_params
#' 
#' @param simplify logical; if \code{TRUE} the output data structure is 
#'   simplified if only one module is specified.
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
#' 
#' @return 
#'  A list of network properties for each module of interest:
#'  \itemize{
#'    \item{\code{'degree'}:}{
#'      The weighted within-module degree: the sum of edge weights for each 
#'      node in the module.
#'    }
#'    \item{\code{'avgWeight'}:}{
#'      The average edge weight within the module.
#'    }
#'  }
#'  If the data used to infer the \code{test} network is provided  
#'  then the following are also returned:
#'  \itemize{
#'    \item{\code{'summary'}:}{
#'      A vector summarising the module across each sample. This is calculated 
#'      as the first eigenvector of the module from a principal component 
#'      analysis.
#'    }
#'    \item{\code{'contribution'}:}{
#'      The \emph{node contribution}: the similarity between each node and the
#'      \emph{module summary profile} (\code{'summary'}).
#'    }
#'    \item{\code{'coherence'}:}{
#'      The proportion of module variance explained by the \code{'summary'}
#'      vector.
#'    }
#'  }
#' 
#' @examples
#' \dontrun{
#' ## Create some example data
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
#' ## Example 1: calculate network properties for a single module
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
#' }
#' 
#' @rdname networkProperties
#' @export
networkProperties <- function(
  data=NULL, correlation, network, moduleAssignments, modules=NULL,
  backgroundLabel="0", discovery=1, test=1, simplify=TRUE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(TRUE, 0, "Validating user input...")
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose=TRUE, tmp.dir)
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
  
  res <- lapply(discovery, function(di) {
    r2 <- lapply(test[[di]], function(ti) {
      r1 <- lapply(modules[[di]], function(mod) {
        
        # Get the row/column indices of the module in the dataset of interest 
        sub <- moduleAssignments[[di]][moduleAssignments[[di]] == mod]
        modInds <- match(names(sub), rownames(network[[ti]]))
        na.inds <- which(is.na(modInds))
        modInds <- na.omit(modInds)
        
        if (length(modInds) == 0) {
          warning(
            'none of the variables composing module "', mod, 
            '" from dataset, "', di, '" are present in dataset "', ti
          )
          return(NULL)
        }
        
        # Get the properties calculated from the underlying data used to infer the
        # network
        datProps <- NULL
        if (!is.null(scaledData[[ti]])) {
          datProps <- dataProps(scaledData[[ti]], modInds)
          # rename for clarity
          names(datProps) <- c("summary", "contribution", "coherence")
          datProps[[2]] <- insert.nas(datProps[[2]], na.inds)
          names(datProps[[1]]) <- rownames(scaledData)
          names(datProps[[2]]) <- names(sub)
        }
        
        # Get the properties calculated from the network.
        netProps <- netProps(network[[ti]], modInds)
        names(netProps) <- c("degree", "avgWeight")
        netProps[[1]] <- insert.nas(netProps[[1]], na.inds)
        names(netProps[[1]]) <- names(sub)
        
        c(datProps, netProps)
      })
      names(r1) <- modules[[di]]
      return(r1)
    })
    names(r2) <- test[[di]]
    return(r2)
  })
  names(res) <- discovery
  
  if (simplify) {
    res <- simplifyList2(res)
  }
  res
}

#' Order nodes and modules within a network.
#' 
#' Order nodes in descending order of \emph{weighted degree} and order 
#' modules by the similarity of their summary vectors.
#' 
#' @inheritParams common_params
#' 
#' @param na.rm logical; If \code{TRUE}, genes present in the \code{discovery} 
#'   dataset but missing from the test dataset are excluded. If \code{FALSE}, 
#'   missing genes are put last in the ordering.
#' @param orderModules logical; if \code{TRUE} modules ordered by clustering
#'   their summary vectors. If \code{FALSE} modules are returned 
#'   in the order provided.
#' @param simplify logical; if \code{FALSE} the returned data structure will be 
#'   a list of vectors, one list element for each module. If \code{TRUE}, the
#'   returned data structure will be a single vector of ordered genes.
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
#'   \link[=plotTopology]{functions}, \link[=sampleOrder]{sample} ordering, and
#'   \link{networkProperties} functions also expect 'bigMatrix' objects. 
#'   Further, 'bigMatrix' objects have a number of benefits, including 
#'   instantaneous load time from any future R session, and parallel access from
#'   mutliple independent R sessions. Methods are provided for 
#'   \link[=bigMatrix-get]{converting to, loading in}, and 
#'   \link[=bigMatrix-out]{writing out} 'bigMatrix' objects.
#' }
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
#'  A vector of variable names in descending order of weighted degree
#'  for each module. 
#'  
#' @examples
#' \dontrun{
#' ## Create some example data
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
#' ## Example 1: get the ordering of samples for a single module
#' nodeOrder(
#'   geA, coexpA, adjA, moduleAssignments, modules="2"
#' )
#' 
#' ## Example 2: get the order of genes of an arbitrary subset
#' ## (the first 10 genes)
#' nodeOrder(
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
#' nodeOrder(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' }
#' 
#' @name nodeOrder
#' @export
nodeOrder <- function(
  data=NULL, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, na.rm=FALSE, orderModules=TRUE, simplify=TRUE
) {
  props <- networkProperties(
    data, correlation, network, moduleAssignments, modules,
    discovery, test, simplify=FALSE
  )
  
  # order modules
  moduleOrder <- 1
  if (length(props) > 1 && orderModules) {
    if (!is.null(data[[test]])) {
      # Create a matrix of module summary vectors to measure the similarity
      seps <- matrix(
        0, ncol=length(props), nrow=length(props[[1]]$summary)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summary
      }
      moduleOrder <- hclust(as.dist(1-cor(seps)))$order
    } else {
      warning(
        "Data used to infer networks not provided, modules will be ordered as", 
        "provided."
      )
      moduleOrder <- seq_along(props)
    }
  } else {
    moduleOrder <- seq_along(props)
  }
  
  # order genes
  res <- foreach(mi = moduleOrder, .combine=c) %do% {
    rr <- names(sort(
      props[[mi]]$degree, decreasing=TRUE, 
      na.last=ifelse(na.rm, NA, TRUE)
    ))
    if (!simplify)
      rr <- list(rr)
    rr
  }
  if (!simplify) {
    names(res) <- names(props)[moduleOrder]
  }
  res
}

#' Order samples within a network.
#' 
#' Get the order of samples within a module based on the module summary vector.
#' 
#' @inheritParams common_params
#'
#' @param na.rm logical; If \code{TRUE} variables present in the 
#'   \code{discovery} dataset but missing from the \code{test} dataset are 
#'   excluded. If \code{FALSE} missing variables are put last in the ordering.
#' @param simplify logical; If \code{TRUE} a vector is returned instead of a 
#'  list when applying the function to one module.
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
#'   \link[=plotTopology]{functions}, \link[=nodeOrder]{nodeOrder} and 
#'   \link{networkProperties} functions also expect 'bigMatrix' objects.
#'   Further, 'bigMatrix' objects have a number of benefits, including 
#'   instantaneous load time from any future R session, and parallel access from
#'   mutliple independent R sessions. Methods are provided for 
#'   \link[=bigMatrix-get]{converting to, loading in}, and 
#'   \link[=bigMatrix-out]{writing out} 'bigMatrix' objects.
#' }
#' 
#' @return
#'  A list of vectors, one per module of interest, each containing the sample 
#'  names sorted in descending order of the module's summary vector. 
#'  
#' @examples
#' \dontrun{
#' ## Create some example data
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
#' ## Example 1: get the ordering of samples for a single module
#' sampleOrder(
#'   geA, coexpA, adjA, moduleAssignments, modules="2"
#' )
#' 
#' ## Example 2: get the order of genes of an arbitrary subset
#' ## (the first 10 genes)
#' sampleOrder(
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
#' # Get the order of samples in the liver tissue for modules 
#' # 3 and 7, which were discovered in the adipose tissue. 
#' sampleOrder(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' }
#' 
#' @name sampleOrder
#' @export
sampleOrder <- function(
  data, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, na.rm=FALSE, simplify=TRUE
) {
  if (is.null(data[[test]]))
    stop("Cannot order samples without 'data'")
  props <- networkProperties(
    data, correlation, network, moduleAssignments, modules,
    discovery, test, simplify=FALSE
  )

  res <- lapply(props, function(mip) {
    # Need to handle cases where no rownames are provided
    if (!is.null(names(mip$summary))) {
      names(sort(
        mip$summary, decreasing=TRUE, 
        na.last=ifelse(na.rm, NA, TRUE)
      ))
    } else {
      order(
        mip$summary, decreasing=TRUE, 
        na.last=ifelse(na.rm, NA, TRUE)
      )
    }
  })
  if (simplify && length(res) == 1) {
    res <- res[[1]]
  }
  res
}
