#' Plot the coexpression of one or more modules
#' 
#' Plot the coexpression of specified modules in either the discovery dataset or
#' another independent test dataset.  If multiple modules are specified, boxes
#' will be drawn around the module assignments.
#' 
#' @param symmetric logical; if \code{TRUE} the coexpression will be plotted as
#'  a symmetric heatmap, if \code{FALSE} it will be plotted as a triangular
#'  heatmap.
#' @param orderBy one of "discovery", "test", or "none". If "discovery" genes
#'   are ordered by intramodular connectivity in the \code{discovery} dataset
#'   (see \code{\link{geneOrder}}). If "test" genes are orderd by intramodular
#'   connectivity in the test dataset. If "none" no ordering is applied.
#' @param plotGeneNames logical; if \code{TRUE}, plot the gene names below the
#'  heatmap.
#' @param plotModuleNames logical; if \code{TRUE}, plot the module names below 
#'   the heatmap. By default, module names are only plotted if multiple
#'   \code{modules} are provided.
#' @param main title for the plot.
#' @param palette a vector of colors to interpolate over when plotting the 
#'   coexpression. The first element should correspond to the color
#'   used when the coexpression between two genes is equal to -1, and the last
#'   element should correspond to the color used when the coexpression between
#'   two genes is equal to 1
#' 
#' @template api_inputs
#' 
#' @examples
#' \dontrun{
#' ## Example 1: Plot coexpression of a single module
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
#' # Plot coexpression of module 2
#' plotCoexpression(
#'   geA, coexpA, adjA, moduleAssignments, modules="2"
#' )
#' 
#' ## Example 2: Plot coexpression of the first 10 genes, 
#' ## without ordering them
#' plotCoexpression(
#'  geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10], orderBy="none"
#' )
#' 
#' ## Example 3: Plot the coexpression of two adipose modules in the liver
#' ## dataset
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
#' # Plot coexpression for adipose modules 3 and 7 in the heart tissue.
#' par(mar=c(7,4,4,2)+0.1)
#' plotCoexpression(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   modules=c("3", "7"), discovery="adipose", test="heart",
#'   main="Heart", plotGeneNames=FALSE
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*") 
#' }
#' 
#' @rdname plotCoexpression
#' @export
plotCoexpression <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, orderBy="discovery", 
  plotGeneNames=TRUE, plotModuleNames, main="Coexpression", 
  palette=coexpression.palette()
) {
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  orderByArgs <- c("discovery", "test", "none")
  orderBy <- orderByArgs[pmatch(orderBy, orderByArgs, nomatch=3)]
  
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
  
  if (orderBy != "none") {
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, test=ifelse(orderBy == "discovery", discovery, test)
    )
  } else {
    geneOrder <- getGenes(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery
    )
  }

 
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  if (symmetric) {
    plotSquareHeatmap(
      coexpression[[test]][geneOrder, geneOrder], palette, vlim=c(-1, 1),
      moduleAssignments[[discovery]][geneOrder]
    )
  } else {
    plotTriangleHeatmap(
      coexpression[[test]][geneOrder, geneOrder], palette, vlim=c(-1, 1),
      moduleAssignments[[discovery]][geneOrder]
    )
  }
  
  # Add axes if specified
  if (plotGeneNames) {
    axis(
      side=1, las=2, at=seq_along(geneOrder), labels=geneOrder, tick=FALSE,
      line=-0.5
    )
  }
  if (plotModuleNames) {
    line <- ifelse(plotGeneNames, 4, -0.5)
    axis(
      side=1, las=1, 
      at=getModuleMidPoints(moduleAssignments[[discovery]][geneOrder]),
      labels=modules, line=line, tick=FALSE
    )
  }
  mtext(main, cex=par("cex.main"), font=2)
}