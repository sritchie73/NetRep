#' Plot a topological feature of network module
#' 
#' Functions to plot individual components of a module's network topology.
#' 
#' @param symmetric logical; if \code{TRUE} the coexpression will be plotted as
#'  a symmetric heatmap, if \code{FALSE} it will be plotted as a triangular
#'  heatmap.
#' @param orderGenesBy one of "discovery", "test", or "none". If "discovery"
#'   genes are ordered by intramodular connectivity in the \code{discovery}
#'   dataset. If "test" genes are orderd by intramodular connectivity in the
#'   \code{test} dataset. If "none" no ordering is applied.
#' @param orderSamplesBy one of "discovery", "test", or "none". If
#'   "discovery"samples are ordered by their summary expression profile in the
#'   \code{discovery} dataset (see \code{\link{sampleOrder}}). If "test" samples
#'   are ordered by their summary expression profile in the \code{test} dataset.
#'   If "none" no ordering is applied.
#' @param orderModules logical; if \code{TRUE} modules ordered by similarity of 
#'   their summary expression profiles. If \code{FALSE} modules are rendered in
#'   the order provided. The default is to order by modules if the gene
#'   expression is provided.
#' @param plotGeneNames logical; if \code{TRUE}, plot the gene names below the
#'  heatmap.
#' @param plotSampleNames logical; if \code{TRUE} the sample names will be 
#'  plotted next to the gene expression heatmap and summary expression profile
#'  plots
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
#'  in \code{plotModuleMembership}, \code{plotConnectivity}, or
#'  \code{plotSummaryExpression}.
#' @param legend logical; if \code{TRUE} legends are drawn for
#'  \code{plotCoexpression}, \code{plotAdjacency}, or \code{plotExpression}.
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
#'  geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10], orderGenesBy="none"
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
plotExpression <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderGenesBy="discovery",
  orderModules, plotGeneNames=TRUE, plotSampleNames=TRUE, 
  plotModuleNames, main="Gene expression", palette=expression.palette(), 
  legend=TRUE
) {
  if (is.null(geneExpression))
    stop("Cannot plot gene expression without gene expression data")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
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
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # Sanity check input for consistency.
  checkSets(
    geneExpression, coexpression, adjacency, moduleAssignments, discovery, test
  )
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot gene expression without gene expression data")
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  
  # Get the summary expression for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Determine gene ordering, then sample ordering.
  if (orderGenesBy == "discovery" && discovery != test) {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    
    propsDisc <- networkProperties(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE
    ) 
    
    moduleOrder <- names(propsDisc)
    if (length(propsDisc) > 1 && orderModules) {
      # Create a matrix of summary expression profiles to measure the similarity
      seps <- matrix(
        0, ncol=length(propsDisc), 
        nrow=length(propsDisc[[1]]$summaryExpression)
      )
      colnames(seps) <- names(propsDisc)
      for (mi in seq_along(propsDisc)) {
        seps[,mi] <- propsDisc[[mi]]$summaryExpression
      }
      moduleOrder <- names(propsDisc)[hclust(as.dist(1-cor(seps)))$order]
    }
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        propsDisc[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  } else if (orderGenesBy == "none") {
    moduleOrder <- names(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    
    # Order modules and samples by the test network
    moduleOrder <- names(props)
    if (length(props) > 1 && orderModules) {
      seps <- matrix(
        0, ncol=length(props), 
        nrow=length(props[[1]]$summaryExpression)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summaryExpression
      }
      moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
    }
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  # If test == discovery, we order samples and modules by discovery, and plot
  # the discovery. So we need to make sure the gene expression exists in the
  # discovery
  if (orderSamplesBy == "discovery" && discovery != test) {
    if (is.null(geneExpression[[discovery]])) {
      stop(
        "Expecting gene expression data for the discovery dataset in order",
        " to sort samples"
      )
    }
    if (!exists("propsDisc")) {
      propsDisc <- networkProperties(
        geneExpression, coexpression, adjacency, moduleAssignments, modules,
        discovery, discovery, FALSE
      )
    }
    sampleOrder <- names(sort(
      propsDisc[[1]]$summaryExpression, decreasing=TRUE
    ))
  } else if (orderSamplesBy == "none") {
    sampleOrder <- rownames(geneExpression[[test]])
  } else {
    sampleOrder <- names(sort(
      props[moduleOrder][[1]]$summaryExpression, decreasing=TRUE
    ))
  }
  
  # Handle missing samples
  if (all(sampleOrder %nin% rownames(geneExpression[[test]]))) {
    stop(
      "No samples from the 'orderSamplesBy' dataset are present in the",
      " 'test' dataset"
    )
  }
  na.pos.y <- which(sampleOrder %nin% rownames(geneExpression[[test]]))
  if (length(na.pos.y) > 0) {
    presentSamples <- sampleOrder[-na.pos.y]
  } else {
    presentSamples <- sampleOrder
  }
  
  # Handle genes not present in the test dataset
  na.pos.x <- which(geneOrder %nin% colnames(coexpression[[test]]))
  if (length(na.pos.x) > 0) {
    presentGenes <- geneOrder[-na.pos.x]
  } else {
    presentGenes <- geneOrder
  }
  
  # Handle expression palette and values
  ge <- geneExpression[[test]][presentSamples, presentGenes]
  range.ge <- range(ge)
  if (all(rr > 0)) {
    palette <- tail(palette, length(palette)/2)
  } else if (all(rr < 0)) {
    palette <- head(palette, length(palette)/2)
  } else {
    range.pal <- c(-max(abs(range.ge)), max(abs(range.ge)))
  }
  plotSquareHeatmap(
    ge, palette, vlim=range.pal,
    moduleAssignments[[discovery]][geneOrder], na.pos.x, na.pos.y
  )
  if (legend) {
    pw <- length(geneOrder) + 1
    ph <- length(sampleOrder) + 1
    addGradientLegend(
      palette, range.pal, range.ge, FALSE, "Expression",
      xlim=c(pw - 0.5 + pw*0.2, pw - 0.5 + pw*0.25), 
      ylim=c(ph/3, ph - 0.5 - ph*0.1)
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
      labels=unique(mas), line=line, tick=FALSE,
      cex.axis=par("cex.lab")
    )
  }
  if (plotSampleNames) {
    axis(
      side=2, las=2, at=rev(seq_along(sampleOrder)), labels=sampleOrder, 
      tick=FALSE, line=-0.5
    )
  }
  mtext(main, cex=par("cex.main"), font=2)
}

#' @rdname plotTopology
#' @export
plotCoexpression <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, orderGenesBy="discovery", orderModules,
  plotGeneNames=TRUE, plotModuleNames, main="Coexpression", 
  palette=coexpression.palette(), legend=TRUE
) {
  if (class(main) != "character")
    stop("'main' must be a characer vector")

  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]

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
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # Get gene ordering
  if (orderGenesBy == "none") {
    geneOrder <- getGenes(moduleAssignments, modules, discovery)
  } else if (orderGenesBy == "discovery") {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    )
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, test, FALSE, orderModules
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
      moduleAssignments[[discovery]][geneOrder], na.pos, na.pos
    )
    if (legend) {
      pd <- length(geneOrder) + 1
      addGradientLegend(
        palette, c(-1,1), c(-1,1), FALSE, "Coexpression",
        xlim=c(pd - 0.5 + pd*0.2, pd - 0.5 + pd*0.25), 
        ylim=c(pd/3, pd - 0.5 - pd*0.1)
      )
    }
  } else {
    plotTriangleHeatmap(
      coexpression[[test]][presentGenes , presentGenes], palette, vlim=c(-1, 1),
      moduleAssignments[[discovery]][geneOrder], na.pos
    )
    if (legend) {
      ph <- length(geneOrder)/2 + 0.5
      pw <- length(geneOrder) + 1
      addGradientLegend(
        palette, c(-1,1), c(-1,1), TRUE, "Coexpression",
        xlim=c(0.5 - pw*0.05, pw*0.22), 
        ylim=c(ph/2 + ph*0.17, ph/2 + ph*0.27)
      )
    }
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
      labels=unique(mas), line=line, tick=FALSE,
      cex.axis=par("cex.lab")
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
        labels=unique(mas), line=line, tick=FALSE,
        cex.axis=par("cex.lab")
      )
    }
  }
  mtext(main, cex=par("cex.main"), font=2)
}

#' @rdname plotTopology
#' @export
plotAdjacency <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, orderGenesBy="discovery", orderModules,
  plotGeneNames=TRUE, plotModuleNames, main="Adjacency", 
  palette=adjacency.palette(), legend=TRUE
) {
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  
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
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # Sanity check input for consistency.
  checkSets(
    geneExpression, coexpression, adjacency, moduleAssignments, discovery, test
  )
  
  # Get gene ordering
  if (orderGenesBy == "none") {
    geneOrder <- getGenes(moduleAssignments, modules, discovery)
  } else if (orderGenesBy == "discovery") {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    )
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, test, FALSE, orderModules
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
      adjacency[[test]][presentGenes, presentGenes], palette, vlim=c(0, 1),
      moduleAssignments[[discovery]][geneOrder], na.pos, na.pos
    )
    if (legend) {
      pd <- length(geneOrder) + 1
      addGradientLegend(
        palette, c(0,1), c(0,1), FALSE, "Adjacency",
        xlim=c(pd - 0.5 + pd*0.2, pd - 0.5 + pd*0.25), 
        ylim=c(pd/3, pd - 0.5 - pd*0.1)
      )
    }
  } else {
    plotTriangleHeatmap(
      adjacency[[test]][presentGenes , presentGenes], palette, vlim=c(0, 1),
      moduleAssignments[[discovery]][geneOrder], na.pos
    )
    if (legend) {
      ph <- length(geneOrder)/2 + 0.5
      pw <- length(geneOrder) + 1
      addGradientLegend(
        palette, c(0,1), c(0,1), TRUE, "Adjacency",
        xlim=c(0.5 - pw*0.05, pw*0.22), 
        ylim=c(ph/2 + ph*0.17, ph/2 + ph*0.27)
      )
    }
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
      labels=unique(mas), line=line, tick=FALSE,
      cex.axis=par("cex.lab")
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
        labels=unique(mas), line=line, tick=FALSE,
        cex.axis=par("cex.lab")
      )
    }
  }
  mtext(main, cex=par("cex.main"), font=2)
}

#' @rdname plotTopology
#' @export
plotModuleMembership <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderGenesBy="discovery", orderModules,
  plotGeneNames=TRUE, plotModuleNames, main="Module Membership", 
  palette=c("#313695", "#a50026"), drawBorder=FALSE
) {
  if (is.null(geneExpression))
    stop("Cannot plot module membership without gene expression data")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  
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
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }

  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot module membership without gene expression data")
  
  # Get the module membership for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )

  # Now we will order the genes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (orderGenesBy == "discovery" && discovery != test)  {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    # Ordering genes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    ) 
  } else if (orderGenesBy == "none") {
    moduleOrder <- seq_along(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- TRUE
    # order modules
    moduleOrder <- 1
    if (length(props) > 1 && orderModules) {
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
      moduleOrder <- seq_along(props)
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
      labels=unique(mas), line=line, tick=FALSE,
      cex.axis=par("cex.lab")
    )
  }
  mtext(main, cex=par("cex.main"), font=2)
}

#' @rdname plotTopology
#' @export
plotConnectivity <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderGenesBy="discovery", orderModules=TRUE,
  plotGeneNames=TRUE, plotModuleNames, main="Normalised Connectivity", 
  palette="#feb24c", drawBorder=FALSE
) {
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  
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
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # Get the connectivity for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Now we will order the genes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (orderGenesBy == "discovery" && discovery != test) {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    # Ordering genes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    ) 
  } else if (orderGenesBy == "none") {
    moduleOrder <- seq_along(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    # order modules
    moduleOrder <- 1
    if (length(props) > 1 && orderModules) {
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
    } else {
      moduleOrder <- seq_along(props)
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
    # Normalise the connectivity by the maximum. The value has no meaning,
    # just the relative sizes and ranks
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
      labels=unique(mas), line=line, tick=FALSE,
      cex.axis=par("cex.lab")
    )
  }
  mtext(main, cex=par("cex.main"), font=2)
}

#' @rdname plotTopology
#' @export
plotSummaryExpression <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderGenesBy="discovery",
  orderModules, plotSampleNames=TRUE, plotModuleNames, 
  main="Summary Expression", palette=c("#762a83", "#1b7837"), drawBorder=FALSE
) {
  if (is.null(geneExpression))
    stop("Cannot plot summary expression without gene expression data")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
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
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot summary expression without gene expression data")
  
  # Get the summary expression for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # We use gene ordering to determine module ordering for consistency with
  # 'plotExpression'. 
  if (orderGenesBy == "discovery" && discovery != test) {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    
    propsDisc <- networkProperties(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE
    ) 
    
    moduleOrder <- names(propsDisc)
    if (length(propsDisc) > 1 && orderModules) {
      # Create a matrix of summary expression profiles to measure the similarity
      seps <- matrix(
        0, ncol=length(propsDisc), 
        nrow=length(propsDisc[[1]]$summaryExpression)
      )
      colnames(seps) <- names(propsDisc)
      for (mi in seq_along(propsDisc)) {
        seps[,mi] <- propsDisc[[mi]]$summaryExpression
      }
      moduleOrder <- names(propsDisc)[hclust(as.dist(1-cor(seps)))$order]
    }
  } else if (orderGenesBy == "none") {
    moduleOrder <- names(props)
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    
    # Order modules and samples by the test network
    moduleOrder <- names(props)
    if (length(props) > 1 && orderModules) {
      seps <- matrix(
        0, ncol=length(props), 
        nrow=length(props[[1]]$summaryExpression)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summaryExpression
      }
      moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
    }
  }
  
  # If test == discovery, we order samples and modules by discovery, and plot
  # the discovery. So we need to make sure the gene expression exists in the
  # discovery
  if (orderSamplesBy == "discovery" && discovery != test) {
    if (is.null(geneExpression[[discovery]])) {
      stop(
        "Expecting gene expression data for the discovery dataset in order",
        " to sort samples"
      )
    }
    if (!exists("propsDisc")) {
      propsDisc <- networkProperties(
        geneExpression, coexpression, adjacency, moduleAssignments, modules,
        discovery, discovery, FALSE
      )
    }
    sampleOrder <- names(sort(
      propsDisc[[1]]$summaryExpression, decreasing=TRUE
    ))
  } else if (orderSamplesBy == "none") {
    sampleOrder <- rownames(geneExpression[[test]])
  } else {
    sampleOrder <- names(sort(
      props[moduleOrder][[1]]$summaryExpression, decreasing=TRUE
    ))
  }
  
  # Handle missing samples
  if (all(sampleOrder %nin% rownames(geneExpression[[test]]))) {
    stop(
      "No samples from the 'orderSamplesBy' dataset are present in the",
      " 'test' dataset"
    )
  }
  na.pos <- which(sampleOrder %nin% rownames(geneExpression[[test]]))
  if (length(na.pos) > 0) {
    presentSamples <- sampleOrder[-na.pos]
  } else {
    presentSamples <- sampleOrder
  }
  
  # now build the Summary Expression matrix
  SEP <- foreach(mi = moduleOrder, .combine=cbind) %do% {
    matrix(
      insert.nas(props[[mi]]$summaryExpression[presentSamples], na.pos),
      ncol=1
    )
  }
  if (plotModuleNames) {
    colnames(SEP) <- moduleOrder
  } else {
    colnames(SEP) <- NULL
  }
  
  # Now build the colors
  cols <- matrix(palette[1], nrow(SEP), ncol(SEP))
  cols[SEP > 0] <- palette[2]
  
  # Plot bar chart
  plotMultiBar(
    SEP, rep(list(range(SEP, na.rm=TRUE)), ncol(SEP)),
    cols=cols, drawBorder=drawBorder, main=main
  )
  
  if (plotSampleNames) {
    axis(
      side=2, tick=FALSE, las=2, at=1:length(sampleOrder)-0.5,
      labels=sampleOrder, line=0
    )
  }
}

#' @rdname plotTopology
#' @export
plotExpressionLegend <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, palette=expression.palette(), main="Expression", 
  horizontal=TRUE
) {
  if (is.null(geneExpression))
    stop("Cannot plot expression legend without gene expression data")
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot summary expression without gene expression data")
  
  modGenes <- getGenes(moduleAssignments, modules, discovery)
  modGenes <- modGenes %sub_in% colnames(geneExpression[[test]])
  if (length(modGenes) == 0)
    stop("None of the module genes are present in the test dataset")
  rg <- range(geneExpression[[test]][,modGenes])
  emptyPlot(c(0,1), c(0,1))
  if (all(rg < 0)) {
    addGradientLegend(head(palette, length(palette)/2), rg, rg, horizontal, main)
  } else if (all(rg > 0)) {
    addGradientLegend(tail(palette, length(palette)/2), rg, rg, horizontal, main)
  } else {
    plim <- c(-max(abs(rg)), max(abs(rg)))
    addGradientLegend(palette, plim, rg, horizontal, main)
  }
}

#' @param horizontal logical; if \code{TRUE} the legend is plotted horizontally.
#' 
#' @rdname plotTopology
#' @export
plotCoexpressionLegend <- function(
  palette=coexpression.palette(), main="Coexpression", horizontal=TRUE
) {
  emptyPlot(c(0,1), c(0,1))
  addGradientLegend(palette, c(-1,1), c(-1,1), horizontal, main)
}

#' @rdname plotTopology
#' @export
plotAdjacencyLegend <- function(
  palette=adjacency.palette(), main="Adjacency", horizontal=TRUE
) {
  emptyPlot(c(0,1), c(0,1))
  addGradientLegend(palette, c(0,1), c(0,1), horizontal, main)
}
