#' Plot a topological feature of network module
#' 
#' 
#' @param symmetric logical; if \code{TRUE} the coexpression will be plotted as
#'  a symmetric heatmap, if \code{FALSE} it will be plotted as a triangular
#'  heatmap.
#' @param geneOrder one of "discovery", "test", "none", or a vector of gene 
#'   names. If "discovery" genes are ordered by intramodular connectivity in the
#'   \code{discovery} dataset (see \code{\link{geneOrder}}). If "test" genes are
#'   orderd by intramodular connectivity in the test dataset. If "none" no 
#'   ordering is applied.
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
#' @param drawBorder logical; if \code{TRUE}, borders are drawn around the bars
#'  in \code{plotModuleMembership}, \code{plotConnectivity}, and 
#'  \code{plotSummaryExpression}.
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
#'  geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10], geneOrder="none"
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
#' @rdname plotTopology
#' @export
plotCoexpression <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, geneOrder="discovery", 
  plotGeneNames=TRUE, plotModuleNames, main="Coexpression", 
  palette=coexpression.palette()
) {
  if (class(main) != "character")
    stop("'main' must be a characer vector")

  geneOrderArgs <- c("discovery", "test", "none")
  geneOrder <- geneOrderArgs[pmatch(geneOrder, geneOrderArgs, nomatch=3)]

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
  
  if (geneOrder == "none") {
    geneOrder <- getGenes(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery
    )
  } else {
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, test=ifelse(geneOrder == "discovery", discovery, test)
    )
  }
 
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  # Handle genes not present in the test dataset
  na.pos <- which(geneOrder %nin% colnames(coexpression[[test]]))
  if (length(na.pos) > 0) {
    presentGenes <- geneOrder[-na.pos]
  } else {
    presentGenes <- geneOrder
  }
  if (symmetric) {
    plotSquareHeatmap(
      coexpression[[test]][presentGenes, presentGenes], palette, vlim=c(-1, 1),
      moduleAssignments[[discovery]][geneOrder], na.pos
    )
  } else {
    plotTriangleHeatmap(
      coexpression[[test]][presentGenes , presentGenes], palette, vlim=c(-1, 1),
      moduleAssignments[[discovery]][geneOrder], na.pos
    )
  }
  
  # Add axes if specified
  mas <- moduleAssignments[[discovery]][geneOrder]
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
      at=getModuleMidPoints(mas),
      labels=modules, line=line, tick=FALSE
    )
  }
  if (symmetric) {
    if (plotGeneNames) {
      axis(
        side=2, las=2, at=rev(seq_along(geneOrder)), labels=geneOrder, tick=FALSE,
        line=-0.5
      )
    }
    if (plotModuleNames) {
      line <- ifelse(plotGeneNames, 4, -0.5)
      axis(
        side=2, las=2, 
        at=length(geneOrder) + 0.5 - getModuleMidPoints(mas),
        labels=modules, line=line, tick=FALSE
      )
    }
  }
  mtext(main, cex=par("cex.main"), font=2)
}

#' @rdname plotTopology
#' @export
plotAdjacency <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, geneOrder="discovery", 
  plotGeneNames=TRUE, plotModuleNames, main="Adjacency", 
  palette=coexpression.palette()
) {
  plotCoexpression(
    geneExpression, adjacency, adjacency, moduleAssignments, modules,
    discovery, test, symmetric, geneOrder, plotGeneNames, plotModuleNames,
    main, palette
  ) 
}

#' @rdname plotTopology
#' @export
plotModuleMembership <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, geneOrder="discovery", 
  plotGeneNames=TRUE, plotModuleNames, main="Module Membership", 
  palette=c("#313695", "#a50026"), drawBorder=FALSE
) {
  if (is.null(geneExpression))
    stop("Cannot plot module membership without gene expression data")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  geneOrderArgs <- c("discovery", "test", "none")
  geneOrder <- geneOrderArgs[pmatch(geneOrder, geneOrderArgs, nomatch=3)]
  
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
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot module membership without gene expression data")
  
  # Get the module membership for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test
  )
  if (length(modules) == 1)
    props <- list(props)
  
  # Now we will order the genes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (geneOrder == "discovery" || (geneOrder == "test" && discovery != test)) {
    # Ordering genes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery
    ) 
  } else if (geneOrder == "none") {
    moduleOrder <- seq_along(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
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
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  # now build the Module Membership vector
  MM <- foreach(mi = seq_along(props), .combine=c) %do% {
    props[[mi]]$moduleMembership
  }
  MM <- MM[geneOrder]
  
  # Plot bar chart
  plotBar(
    MM, c(-1,1), moduleAssignments[[discovery]][geneOrder],
    ifelse(MM > 0, palette[2], palette[1]), drawBorder=drawBorder
  )
  
  # Add axes if specified
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  mas <- moduleAssignments[[discovery]][geneOrder]
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
      at=getModuleMidPoints(mas),
      labels=modules, line=line, tick=FALSE
    )
  }
  mtext(main, cex=par("cex.main"), font=2)
}

#' @rdname plotTopology
#' @export
plotConnectivity <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, geneOrder="discovery", 
  plotGeneNames=TRUE, plotModuleNames, main="Normalised Connectivity", 
  palette="#feb24c", drawBorder=FALSE
) {
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  geneOrderArgs <- c("discovery", "test", "none")
  geneOrder <- geneOrderArgs[pmatch(geneOrder, geneOrderArgs, nomatch=3)]
  
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
  
  # Get the module membership for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test
  )
  if (length(modules) == 1)
    props <- list(props)
  
  # Now we will order the genes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (geneOrder == "discovery" || (geneOrder == "test" && discovery != test)) {
    # Ordering genes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery
    ) 
  } else if (geneOrder == "none") {
    moduleOrder <- seq_along(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
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
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  # now build the (Normalised) Intramodular Connectivity vector
  kIM <- foreach(mi = seq_along(props), .combine=c) %do% {
    props[[mi]]$connectivity/max(na.omit(props[[mi]]$connectivity))
  }
  kIM <- kIM[geneOrder]
  
  # Plot bar chart
  plotBar(
    kIM, c(0,1), moduleAssignments[[discovery]][geneOrder],
    palette, drawBorder=drawBorder
  )
  
  # Add axes if specified
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  mas <- moduleAssignments[[discovery]][geneOrder]
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
      at=getModuleMidPoints(mas),
      labels=modules, line=line, tick=FALSE
    )
  }
  mtext(main, cex=par("cex.main"), font=2)
}
