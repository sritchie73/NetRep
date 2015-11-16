#' Plot a topological feature of network module
#' 
#' Functions for plotting the topology of a network module.
#' 
#' @inheritParams common_params
#' @inheritParams common_params2
#'
#' @param symmetric logical; controls whether the correlation and network 
#'  heatmaps are drawn as symmetric (square) heatmaps or asymettric triangle 
#'  heatmaps. If symmetric, then the node and module names will also be rendered
#'  on the left axis.
#' @param orderNodesBy one of "discovery", "test" or "none". Controls how nodes
#' are ordered on the plot (see details).
#' @param orderSamplesBy one of "discovery", "test" or "none". Controls how 
#' samples are ordered on the plot (see details).
#' @param orderModules logical; if \code{TRUE} modules will be ordered by 
#'  similarity (see details). The default is \code{TRUE} if the data for the 
#'  \code{test} network is provided.
#' @param plotNodeNames logical; controls whether the node names are 
#'  rendered on the bottom axis.
#' @param plotSampleNames logical; controls whether the sample names are 
#'  rendered on the left axis.
#' @param plotModuleNames logical; controls whether module names are rendered.
#'  The default is for module names to be rendered when multiple \code{modules} 
#'  are drawn.
#' @param palette a vector of colors to use for each plot (see details).
#' @param drawBorders logical; if \code{TRUE}, borders are drawn around the bars
#'  in \code{plotModuleMembership}, \code{plotConnectivity}, and
#'  \code{plotModuleSummary}.
#' @param border.width line width for borders.
#' @param plotLegend logical; controls whether a legend is drawn when using
#'  \code{plotCorrelation}, \code{plotNetwork}, or \code{plotData}.
#' @param gaxt.line the number of lines into the bottom margin at which the node
#'  names will be drawn.
#' @param saxt.line the number of lines into the left margin at which the sample
#'  names will be drawn.
#' @param maxt.line the number of lines into the bottom margin at which the 
#'  module names will be drawn.
#' @param legend.main title for the legend.
#' @param main title for each plot.
#' @param legend.position the distance from the plot to start the legend, as a
#'  proportion of the plot width.
#' @param legend.tick.size size of the ticks on the axis legend.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{legend.tick.size}.
#' @param cex.axis relative size of the node and sample names.
#' @param cex.lab relative size of the module names and legend titles.
#' @param cex.main relative size of the plot titles.
#' @param horizontal logical; controls whether the legend is rendered 
#'  horizontally or vertically when using \code{plotDataLegend},
#'  \code{plotCorrelationLegend} or \code{plotNetworkLegend}.
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
#'   \link[=sampleOrder]{sample} ordering functions  also expect 'bigMatrix'
#'   objects. Further, 'bigMatrix' objects have a number of benefits, including 
#'   instantaneous load time from any future R session, and parallel access from
#'   mutliple independent R sessions. Methods are provided for 
#'   \link[=bigMatrix-get]{converting to, loading in}, and 
#'   \link[=bigMatrix-out]{writing out} 'bigMatrix' objects.
#' }
#' \subsection{Node, sample, and module ordering:}{
#'   By default, nodes are ordered in decreasing order of within-module
#'   connectivity in the \code{discovery} dataset (see \code{\link{nodeOrder}}). 
#'   This facilitates the visual comparison of modules across datasets, as the 
#'   node ordering will be preserved. Missing nodes are colored in grey. If
#'   \code{orderNodesBy} is "test" nodes will instead be ordered by 
#'   within-module connectivity in the \code{test} dataset. If "none" nodes are 
#'   drawn in the order they are provided in the drawn dataset.
#'   
#'   When multiple modules are specified, modules are ordered by the similarity
#'   of their summary vectors in the drawn dataset. To disable this behaviour, 
#'   set \code{orderModules} to \code{FALSE}.
#'   
#'   Sample ordering only applies to \code{plotData} and 
#'   \code{plotModuleSummary}. By default, samples are ordered in descending
#'   order of the module summary vector in the drawn dataset for the left-most 
#'   module appearing on the plot (see \code{\link{sampleOrder}}.
#' }
#' \subsection{Normalised connectivity:}{
#'   The within-module connectivity is normalised by the maximum connectivity in
#'   any given module when rendered on the bar plot. This facilitates visual 
#'   comparison of multiple modules with differing sizes or densities.
#' }
#' \subsection{Customising plot layout:}{
#'   Although reasonable default values for most parameters have been provided,
#'   the rendering of axes and titles may need adjusting depending on the size
#'   of the plot window. The parameters \code{gaxt.line}, \code{saxt.line}, 
#'   \code{maxt.line}, and \code{laxt.line} control the distance from each plot
#'   window that the node labels, sample labels, module labels, and legend 
#'   labels are rendered. 
#'   
#'   \code{legend.tick.size} controls the length of the 
#'   axis ticks on each of the legends relative to the correlation, network,
#'   and data plot windows. 
#'   
#'   \code{legend.position} controls the horizontal offset of the legend 
#'   relative to the plot. For the triangle heatmaps, (\code{symmetric=FALSE} in
#'   \code{plotCorrelation} and \code{plotNetwork}) this controls how far 
#'   left of the plot the legend starts as a proportion of the plot width. For
#'   the square heatmaps (\code{plotData}, and \code{symmetric=TRUE} in
#'   \code{plotCorrelation} and \code{plotNetwork}) this controls how far 
#'   right of the plot the legend starts as a proportion of the plot width.
#'   
#'   \code{cex.main} controls the relative text size of the plot title
#'   (specified by the \code{main} argument). \code{cex.axis} controls the
#'   relative text size of the node and sample labels. \code{cex.lab} controls
#'   the relative text size of the bar plot axis labels, module labels, and the
#'   legend titles.
#'   
#'   The rendering of node, sample, and module names can be disabled by setting
#'   \code{plotNodeNames}, \code{plotSampleNames}, and \code{plotModuleNames} to
#'   \code{FALSE}, and the rendering of the legend can be disabled by setting
#'   \code{plotLegend} to \code{FALSE}
#'   
#'   The \code{drawBorders} argument controls whether borders are drawn around
#'   the bars in \code{plotConnectivity}, \code{plotModuleMembership}, and 
#'   \code{plotModuleSummary}.
#' }
#' \subsection{Customising the color palette:}{
#'   \code{plotCorrelation} and \code{plotCorrelationLegend} expect the
#'   \code{palette} argument to be a vector of colors to interpolate over when
#'   plotting the correlation. They expect the first element of the 
#'   \code{palette} vector to be the color used for correlation values of -1,
#'   and the last element of the \code{palette} vector to be the color used for 
#'   correlation values of 1.
#'   
#'   \code{plotNetwork} and \code{plotNetworkLegend} expect the
#'   \code{palette} argument to be a vector of colors to interpolate over when
#'   plotting the edge weights. They expect the first element of the 
#'   \code{palette} vector to be the color used for edge weights of 0,
#'   and the last element of the \code{palette} vector to be the color used for 
#'   correlation values of 1.
#'   
#'   \code{plotConnectivity} expects \code{palette} to be a single color, a 
#'   vector of colors, one for each node, or a vector of colors to be repeated.
#'   
#'   \code{plotModuleMembership} expects \code{palette} to be a vector 
#'   containing two colors, the first to be used for nodes with negative module
#'   membership values, and the second to be used for nodes with positive module
#'   membership values. 
#'   
#'   \code{plotData} and \code{plotDataLegend} expect the \code{palette}
#'   argument to be a vector of colors to interpolate over when plotting the
#'   'data.' In order to accomodate data matrices with different ranges,
#'   these functions expect the palette to be a diverging set of colors with a
#'   centre value of 0 (e.g. white). The colors drawn are balanced around 0: 
#'   i.e. positive and negative values of data will have the same intensity on a 
#'   diverging color palette regardless of the actual range of the data.
#'   
#'   \code{plotModuleSummary} expects \code{palette} to be a vector 
#'   containing two colors, the first to be used for nodes with a negative 
#'   module summary measurement, and the second to be used for genes with a 
#'   positive module summary measurment, regardless of whether the data matrix
#'   is centred around 0.
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
#' @seealso
#' \code{\link{plotModule}} for a combined plot showing all topological 
#' properties for a network module.
#' 
#' @examples
#' \dontrun{
#' ## Create some example data
#' geA <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geA) <- paste0("Gene_", 1:100)
#' rownames(geA) <- paste0("CohortA_", 1:50)
#' coexpA <- cor(geA) # correlation
#' adjA <- abs(coexpA)^5 # network
#' moduleAssignments <- sample(1:7, size=100, replace=TRUE)
#' names(moduleAssignments) <- paste0("Gene_", 1:100)
#' 
#' # Create bigMatrix objects for each matrix.
#' geA <- as.bigMatrix(geA, "geA_bm")
#' coexpA <- as.bigMatrix(coexpA, "coexpA_bm")
#' adjA <- as.bigMatrix(adjA, "adjA_bm")
#' 
#' ## Example 1: Plot Module 2 in cohort A.
#' plotData(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotCorrelation(geA, coexpA, adjA, moduleAssignments, modules="2")
#' # alternatively as a square heatmap
#' plotCorrelation(
#'  geA, coexpA, adjA, moduleAssignments, modules="2", symmetric=TRUE
#' )
#' plotNetwork(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotConnectivity(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotModuleMembership(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotModuleSummary(geA, coexpA, adjA, moduleAssignments, modules="2")
#' 
#' ## Example 2: Plot an arbitrary set of genes in cohort A
#' plotData(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotCorrelation(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' # alternatively as a square heatmap
#' plotCorrelation(
#'  geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10], symmetric=TRUE
#' )
#' plotNetwork(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotConnectivity(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotModuleMembership(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotModuleSummary(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' 
#' ## Example 3: Plot the topology of two adipose tissue modules in the liver
#' ## tissue data 
#' 
#' geAdipose <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geAdipose) <- paste0("Gene_", 1:100)
#' rownames(geAdipose) <- paste0("Sample_", 1:50)
#' coexpAdipose <- cor(geAdipose) # correlation
#' adjAdipose <- abs(coexpAdipose)^5 # network
#' adiposeModules <- sample(0:7, size=100, replace=TRUE)
#' names(adiposeModules) <- paste0("Gene_", 1:100)
#' 
#' geLiver <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geLiver) <- paste0("Gene_", 1:100)
#' rownames(geLiver) <- paste0("Sample_", 1:50)
#' coexpLiver <- cor(geLiver) # correlation
#' adjLiver <- abs(coexpLiver)^6 # network
#' liverModules <- sample(0:12, size=100, replace=TRUE)
#' names(liverModules) <- paste0("Gene_", 1:100)
#'
#' geHeart <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geHeart) <- paste0("Gene_", 1:100)
#' rownames(geHeart) <- paste0("Sample_", 1:50)
#' coexpHeart <- cor(geHeart) # correlation
#' adjHeart <- abs(coexpHeart)^4 # network
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
#' correlation <- list(
#'   adipose=as.bigMatrix(coexpAdipose, "coexpAdipose_bm"),
#'   liver=as.bigMatrix(coexpLiver, "coexpLiver_bm"),  
#'   heart=as.bigMatrix(coexpHeart, "coexpHeart_bm") 
#' )
#' network <- list(
#'   adipose=as.bigMatrix(adjAdipose, "adjAdipose_bm"),
#'   liver=as.bigMatrix(adjLiver, "adjLiver_bm"),  
#'   heart=as.bigMatrix(adjHeart, "adjHeart_bm") 
#' )
#' moduleAssignments <- list(
#'   adipose=adiposeModules, liver=liverModules, heart=heartModules
#' )
#' 
#' # Show the plot
#' plotData(
#'  geneExpression, correlation, network, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotCorrelation(
#'  geneExpression, correlation, network, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' # alternatively as a square heatmap
#' plotCorrelation(
#'  geneExpression, correlation, network, moduleAssignments, 
#'  modules=c("3", "7"), discovery="adipose", test="liver", symmetric=TRUE
#' )
#' plotNetwork(
#'  geneExpression, correlation, network, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotConnectivity(
#'  geneExpression, correlation, network, moduleAssignments, 
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotModuleMembership(
#'  geneExpression, correlation, network, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotModuleSummary(
#'  geneExpression, correlation, network, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' }
#' 
#' @name plotTopology
NULL

#' \code{plotData}: Plot a heatmap of the data matrix for one or more
#' network modules.
#' 
#' @rdname plotTopology
#' @export
plotData <- function(
  data, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderNodesBy="discovery",
  orderModules, plotNodeNames=TRUE, plotSampleNames=TRUE, plotModuleNames,
  main="", palette=data.palette(), border.width=2, plotLegend=TRUE, 
  legend.main="data", gaxt.line=-0.5, saxt.line=-0.5, maxt.line=3, 
  legend.position=0.15, legend.tick.size=0.03, laxt.line=3, cex.axis=0.8, 
  cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (is.null(data))
    stop("Cannot plot data matrix without 'data'")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  data <- unifyDS(dynamicMatLoad(data))
  correlation <- unifyDS(dynamicMatLoad(correlation))
  network <- unifyDS(dynamicMatLoad(network))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  data <- formatDataList(
    data, length(correlation), names(correlation)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(correlation), names(correlation),
    ncol(correlation[[discovery]]), colnames(correlation[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Sanity check input for consistency.
  checkSets(
    data, correlation, network, moduleAssignments, discovery, test
  )
  
  if (is.null(data[[test]]))
    stop("Cannot plot data matrix without 'data'")
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and samples in the 'test' dataset by the dataset 
  # specified in 'orderNodesBy' and 'orderSamplesBy'.
  #-----------------------------------------------------------------------------
  props <- networkProperties(
    data, correlation, network, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Determine node ordering, then sample ordering.
  if (orderNodesBy == "discovery" && discovery != test) {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[discovery]]), FALSE, TRUE)
    
    propsDisc <- networkProperties(
      data, correlation, network, moduleAssignments, modules,
      discovery, discovery, FALSE
    ) 
    
    moduleOrder <- names(propsDisc)
    if (length(propsDisc) > 1 && orderModules) {
      # Create a matrix of module summary vectors to measure the similarity
      seps <- matrix(
        0, ncol=length(propsDisc), 
        nrow=length(propsDisc[[1]]$moduleSummary)
      )
      colnames(seps) <- names(propsDisc)
      for (mi in seq_along(propsDisc)) {
        seps[,mi] <- propsDisc[[mi]]$moduleSummary
      }
      moduleOrder <- names(propsDisc)[hclust(as.dist(1-cor(seps)))$order]
    }
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        propsDisc[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  } else if (orderNodesBy == "none") {
    moduleOrder <- names(props)
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[test]]), FALSE, TRUE)
    
    # Order modules and samples by the test network
    moduleOrder <- names(props)
    if (length(props) > 1 && orderModules) {
      seps <- matrix(
        0, ncol=length(props), 
        nrow=length(props[[1]]$moduleSummary)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$moduleSummary
      }
      moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
    }
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  # If test == discovery, we order samples and modules by discovery, and plot
  # the discovery. So we need to make sure the data matrix exists in the
  # discovery
  if (orderSamplesBy == "discovery" && discovery != test) {
    if (is.null(data[[discovery]])) {
      stop(
        "Expecting data matrix for the discovery dataset in order",
        " to sort samples"
      )
    }
    if (!exists("propsDisc")) {
      propsDisc <- networkProperties(
        data, correlation, network, moduleAssignments, modules,
        discovery, discovery, FALSE
      )
    }
    sampleOrder <- names(sort(
      propsDisc[[1]]$moduleSummary, decreasing=TRUE
    ))
  } else if (orderSamplesBy == "none") {
    sampleOrder <- rownames(data[[test]])
  } else {
    sampleOrder <- names(sort(
      props[moduleOrder][[1]]$moduleSummary, decreasing=TRUE
    ))
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  if (all(sampleOrder %nin% rownames(data[[test]]))) {
    stop(
      "No samples from the 'orderSamplesBy' dataset are present in the",
      " 'test' dataset"
    )
  }
  na.pos.y <- which(sampleOrder %nin% rownames(data[[test]]))
  if (length(na.pos.y) > 0) {
    presentSamples <- sampleOrder[-na.pos.y]
  } else {
    presentSamples <- sampleOrder
  }
  
  # Handle nodes not present in the test dataset
  na.pos.x <- which(nodeOrder %nin% colnames(correlation[[test]]))
  if (length(na.pos.x) > 0) {
    presentVars <- nodeOrder[-na.pos.x]
  } else {
    presentVars <- nodeOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the data matrix
  #-----------------------------------------------------------------------------
  # First we need to set up the color palette for the data, which 
  # includes the fact that the range may be unbalanced around 0.
  ge <- data[[test]][presentSamples, presentVars]
  range.ge <- range(ge)
  if (all(range.ge > 0)) {
    palette <- tail(palette, length(palette)/2)
  } else if (all(range.ge < 0)) {
    palette <- head(palette, length(palette)/2)
  } else {
    range.pal <- c(-max(abs(range.ge)), max(abs(range.ge)))
  }
  xaxt <- NULL
  if (plotNodeNames)
    xaxt <- nodeOrder
  yaxt <- NULL
  if (plotSampleNames)
    yaxt <- sampleOrder
  plotSquareHeatmap(
    ge, palette, vlim=range.pal, legend.lim=range.ge,
    moduleAssignments[[discovery]][nodeOrder], na.pos.x, na.pos.y, 
    xaxt=xaxt, yaxt=yaxt, plotLegend=plotLegend, main=main,
    legend.main=legend.main, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, yaxt.line=saxt.line, legend.tick.size=legend.tick.size,
    laxt.line=laxt.line, legend.line=legend.position, maxt.line=maxt.line,
    border.width=border.width
  )
}

#' \code{plotCorrelation}: Plot a heatmap of the correlation structure for one 
#' or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotCorrelation <- function(
  data=NULL, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, orderNodesBy="discovery", orderModules,
  plotNodeNames=TRUE, plotModuleNames, main="", palette=correlation.palette(), 
  border.width=2, plotLegend=TRUE, legend.main="correlation",
  gaxt.line=-0.5, maxt.line=3, legend.position, legend.tick.size=0.03, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (class(main) != "character")
    stop("'main' must be a characer vector")

  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]

  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  data <- unifyDS(dynamicMatLoad(data))
  correlation <- unifyDS(dynamicMatLoad(correlation))
  network <- unifyDS(dynamicMatLoad(network))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  data <- formatDataList(
    data, length(correlation), names(correlation)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(correlation), names(correlation),
    ncol(correlation[[discovery]]), colnames(correlation[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes in the 'test' dataset by the dataset specified in 
  # 'orderNodesBy'.
  #-----------------------------------------------------------------------------
  if (orderNodesBy == "none") {
    nodeOrder <- getModuleVarsUnsorted(moduleAssignments, modules, discovery)
  } else if (orderNodesBy == "discovery") {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[discovery]]), FALSE, TRUE)
    nodeOrder <- nodeOrder(
      data, correlation, network, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    )
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[test]]), FALSE, TRUE)
    nodeOrder <- nodeOrder(
      data, correlation, network, moduleAssignments, modules,
      discovery, test, FALSE, orderModules
    )
  }
 
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Identify nodes from the 'discovery' dataset not present in the 'test' 
  # dataset.
  #-----------------------------------------------------------------------------
  na.pos <- which(nodeOrder %nin% colnames(correlation[[test]]))
  if (length(na.pos) > 0) {
    presentVars <- nodeOrder[-na.pos]
  } else {
    presentVars <- nodeOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the correlation structure
  #-----------------------------------------------------------------------------
  gaxt <- NULL
  if (plotNodeNames)
    gaxt <- nodeOrder
  if (symmetric) {
    if (missing(legend.position))
      legend.position <- 0.2
    plotSquareHeatmap(
      correlation[[test]][presentVars, presentVars], palette, c(-1, 1), 
      moduleAssignments[[discovery]][nodeOrder], na.pos, na.pos, 
      xaxt=gaxt, yaxt=gaxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, yaxt.line=gaxt.line, border.width=border.width,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line
    )
  } else {
    if (missing(legend.position))
      legend.position <- 0.1
    plotTriangleHeatmap(
      correlation[[test]][presentVars , presentVars], palette, c(-1, 1),
      moduleAssignments[[discovery]][nodeOrder], na.pos, xaxt=gaxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=gaxt.line,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line, 
      border.width=border.width
    )
  }
}

#' \code{plotNetwork}: Plot a heatmap of the edge weights for one or more
#' network modules.
#' 
#' @rdname plotTopology
#' @export
plotNetwork <- function(
  data=NULL, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, orderNodesBy="discovery", orderModules,
  plotNodeNames=TRUE, plotModuleNames, main="", palette=network.palette(), 
  border.width=2, plotLegend=TRUE, legend.main="network",
  gaxt.line=-0.5, maxt.line=3, legend.position, legend.tick.size=0.03, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  data <- unifyDS(dynamicMatLoad(data))
  correlation <- unifyDS(dynamicMatLoad(correlation))
  network <- unifyDS(dynamicMatLoad(network))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  data <- formatDataList(
    data, length(correlation), names(correlation)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(correlation), names(correlation),
    ncol(correlation[[discovery]]), colnames(correlation[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes in the 'test' dataset by the dataset specified in 
  # 'orderNodesBy'.
  #-----------------------------------------------------------------------------
  if (orderNodesBy == "none") {
    nodeOrder <- getModuleVarsUnsorted(moduleAssignments, modules, discovery)
  } else if (orderNodesBy == "discovery") {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[discovery]]), FALSE, TRUE)
    nodeOrder <- nodeOrder(
      data, correlation, network, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    )
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[test]]), FALSE, TRUE)
    nodeOrder <- nodeOrder(
      data, correlation, network, moduleAssignments, modules,
      discovery, test, FALSE, orderModules
    )
  }
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Identify nodes from the 'discovery' dataset not present in the 'test' 
  # dataset.
  #-----------------------------------------------------------------------------
  na.pos <- which(nodeOrder %nin% colnames(correlation[[test]]))
  if (length(na.pos) > 0) {
    presentVars <- nodeOrder[-na.pos]
  } else {
    presentVars <- nodeOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the network edge weights
  #-----------------------------------------------------------------------------
  gaxt <- NULL
  if (plotNodeNames)
    gaxt <- nodeOrder
  if (symmetric) {
    if (missing(legend.position))
      legend.position <- 0.2
    plotSquareHeatmap(
      network[[test]][presentVars, presentVars], palette, c(0, 1), 
      moduleAssignments[[discovery]][nodeOrder], na.pos, na.pos, 
      xaxt=gaxt, yaxt=gaxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, yaxt.line=gaxt.line, 
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line,
      border.width=border.width
    )
  } else {
    if (missing(legend.position))
      legend.position <- 0.1
    plotTriangleHeatmap(
      network[[test]][presentVars , presentVars], palette, c(0, 1),
      moduleAssignments[[discovery]][nodeOrder], na.pos, xaxt=gaxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=gaxt.line,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line, 
      border.width=border.width
    )
  }
}

#' \code{plotModuleMembership}: Plot a bar chart of the module membership for
#' one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotModuleMembership <- function(
  data=NULL, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, orderNodesBy="discovery", orderModules,
  plotNodeNames=TRUE, plotModuleNames, main="", border.width=2, 
  palette=c("#313695", "#a50026"), drawBorders=FALSE, gaxt.line=-0.5, 
  maxt.line=3, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (is.null(data))
    stop("Cannot plot module membership without the corresponding 'data'")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  data <- unifyDS(dynamicMatLoad(data))
  correlation <- unifyDS(dynamicMatLoad(correlation))
  network <- unifyDS(dynamicMatLoad(network))
  
  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(correlation), names(correlation),
    ncol(correlation[[discovery]]), colnames(correlation[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  if (is.null(data[[test]]))
    stop("Cannot plot module membership without the corresponding 'data'")
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes in the 'test' dataset by the dataset specified in 
  # 'orderNodesBy'.
  #-----------------------------------------------------------------------------
  # Get the module membership for each module in the test network.
  props <- networkProperties(
    data, correlation, network, moduleAssignments, modules, 
    discovery, test, FALSE
  )

  # Now we will order the nodes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (orderNodesBy == "discovery" && discovery != test)  {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[discovery]]), FALSE, TRUE)
    # Ordering nodes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    nodeOrder <- nodeOrder(
      data, correlation, network, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    ) 
  } else if (orderNodesBy == "none") {
    moduleOrder <- seq_along(props)
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- TRUE
    # order modules
    moduleOrder <- 1
    if (length(props) > 1 && orderModules) {
      # Create a matrix of module summary vectors to measure the similarity
      seps <- matrix(
        0, ncol=length(props), nrow=length(props[[1]]$moduleSummary)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$moduleSummary
      }
      moduleOrder <- hclust(as.dist(1-cor(seps)))$order
    } else {
      moduleOrder <- seq_along(props)
    }
    
    # order nodes
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Module Membership
  #-----------------------------------------------------------------------------
  # now build the Module Membership vector
  MM <- foreach(mi = seq_along(props), .combine=c) %do% {
    props[[mi]]$moduleMembership
  }
  MM <- MM[nodeOrder]
  
  # Plot bar chart
  plotBar(
    MM, c(-1,1), moduleAssignments[[discovery]][nodeOrder],
    ifelse(MM > 0, palette[2], palette[1]), drawBorders=drawBorders,
    xaxt=plotNodeNames, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, maxt.line=maxt.line, main=main,
    ylab="Module membership", border.width=border.width
  )
}

#' \code{plotConnectivity:} Plot a bar chart of the normalised intramodular 
#' connectivity (see details) for one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotConnectivity <- function(
  data=NULL, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, orderNodesBy="discovery", orderModules=TRUE,
  plotNodeNames=TRUE, plotModuleNames, main="", palette="#feb24c", 
  border.width=2,  drawBorders=FALSE, gaxt.line=-0.5, maxt.line=3, 
  cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  data <- unifyDS(dynamicMatLoad(data))
  correlation <- unifyDS(dynamicMatLoad(correlation))
  network <- unifyDS(dynamicMatLoad(network))
  
  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(correlation), names(correlation),
    ncol(correlation[[discovery]]), colnames(correlation[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes in the 'test' dataset by the dataset specified in 
  # 'orderNodesBy'.
  #-----------------------------------------------------------------------------
  # Get the module membership for each module in the test network.
  props <- networkProperties(
    data, correlation, network, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Now we will order the nodes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (orderNodesBy == "discovery" && discovery != test)  {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[discovery]]), FALSE, TRUE)
    # Ordering nodes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    nodeOrder <- nodeOrder(
      data, correlation, network, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    ) 
  } else if (orderNodesBy == "none") {
    moduleOrder <- seq_along(props)
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- TRUE
    # order modules
    moduleOrder <- 1
    if (length(props) > 1 && orderModules && !is.null(data[[test]])) {
      # Create a matrix of module summary vectors to measure the similarity
      seps <- matrix(
        0, ncol=length(props), nrow=length(props[[1]]$moduleSummary)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$moduleSummary
      }
      moduleOrder <- hclust(as.dist(1-cor(seps)))$order
    } else {
      moduleOrder <- seq_along(props)
    }
    
    # order nodes
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Connectivity
  #-----------------------------------------------------------------------------
  # now build the (Normalised) Intramodular Connectivity vector
  kIM <- foreach(mi = seq_along(props), .combine=c) %do% {
    # Normalise the connectivity by the maximum. The value has no meaning,
    # just the relative sizes and ranks
    props[[mi]]$connectivity/max(na.omit(props[[mi]]$connectivity))
  }
  kIM <- kIM[nodeOrder]
  
  # Plot bar chart
  plotBar(
    kIM, c(0,1), moduleAssignments[[discovery]][nodeOrder],
    palette, drawBorders=drawBorders,
    xaxt=plotNodeNames, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, maxt.line=maxt.line, main=main,
    ylab="Normalised connectivity", border.width=border.width
  )
}

#' \code{plotModuleSummary}: Plot bar charts of the module summary vectors of 
#' one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotModuleSummary <- function(
  data, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderNodesBy="discovery",
  orderModules, plotSampleNames=TRUE, plotModuleNames, main="", 
  palette=c("#762a83", "#1b7837"), border.width=2, drawBorders=FALSE, 
  saxt.line=-0.5, maxt.line=0, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (is.null(data))
    stop("Cannot plot data matrix without 'data'")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  data <- unifyDS(dynamicMatLoad(data))
  correlation <- unifyDS(dynamicMatLoad(correlation))
  network <- unifyDS(dynamicMatLoad(network))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  data <- formatDataList(
    data, length(correlation), names(correlation)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(correlation), names(correlation),
    ncol(correlation[[discovery]]), colnames(correlation[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Sanity check input for consistency.
  checkSets(
    data, correlation, network, moduleAssignments, discovery, test
  )
  
  if (is.null(data[[test]]))
    stop("Cannot plot data matrix without 'data'")
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and samples in the 'test' dataset by the dataset 
  # specified in 'orderNodesBy' and 'orderSamplesBy'.
  #-----------------------------------------------------------------------------
  props <- networkProperties(
    data, correlation, network, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Determine node ordering, then sample ordering.
  if (orderNodesBy == "discovery" && discovery != test) {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[discovery]]), FALSE, TRUE)
    
    propsDisc <- networkProperties(
      data, correlation, network, moduleAssignments, modules,
      discovery, discovery, FALSE
    ) 
    
    moduleOrder <- names(propsDisc)
    if (length(propsDisc) > 1 && orderModules) {
      # Create a matrix of module summary vectors to measure the similarity
      seps <- matrix(
        0, ncol=length(propsDisc), 
        nrow=length(propsDisc[[1]]$moduleSummary)
      )
      colnames(seps) <- names(propsDisc)
      for (mi in seq_along(propsDisc)) {
        seps[,mi] <- propsDisc[[mi]]$moduleSummary
      }
      moduleOrder <- names(propsDisc)[hclust(as.dist(1-cor(seps)))$order]
    }
  } else if (orderNodesBy == "none") {
    moduleOrder <- names(props)
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(data[[test]]), FALSE, TRUE)
    
    # Order modules and samples by the test network
    moduleOrder <- names(props)
    if (length(props) > 1 && orderModules) {
      seps <- matrix(
        0, ncol=length(props), 
        nrow=length(props[[1]]$moduleSummary)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$moduleSummary
      }
      moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
    }
  }
  
  # If test == discovery, we order samples and modules by discovery, and plot
  # the discovery. So we need to make sure the data exists in the
  # discovery
  if (orderSamplesBy == "discovery" && discovery != test) {
    if (is.null(data[[discovery]])) {
      stop(
        "Expecting data for the discovery dataset in order",
        " to sort samples"
      )
    }
    if (!exists("propsDisc")) {
      propsDisc <- networkProperties(
        data, correlation, network, moduleAssignments, modules,
        discovery, discovery, FALSE
      )
    }
    sampleOrder <- names(sort(
      propsDisc[[1]]$moduleSummary, decreasing=TRUE
    ))
  } else if (orderSamplesBy == "none") {
    sampleOrder <- rownames(data[[test]])
  } else {
    sampleOrder <- names(sort(
      props[moduleOrder][[1]]$moduleSummary, decreasing=TRUE
    ))
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  if (all(sampleOrder %nin% rownames(data[[test]]))) {
    stop(
      "No samples from the 'orderSamplesBy' dataset are present in the",
      " 'test' dataset"
    )
  }
  na.pos <- which(sampleOrder %nin% rownames(data[[test]]))
  if (length(na.pos) > 0) {
    presentSamples <- sampleOrder[-na.pos]
  } else {
    presentSamples <- sampleOrder
  }
  
  
  #-----------------------------------------------------------------------------
  # Plot the module summary vectors 
  #-----------------------------------------------------------------------------
  SEP <- foreach(mi = moduleOrder, .combine=cbind) %do% {
    matrix(
      insert.nas(props[[mi]]$moduleSummary[presentSamples], na.pos),
      ncol=1
    )
  }
  colnames(SEP) <- moduleOrder
  rownames(SEP) <- sampleOrder
  
  # Now build the colors
  cols <- matrix(palette[1], nrow(SEP), ncol(SEP))
  cols[SEP > 0] <- palette[2]
  
  # Plot bar chart
  plotMultiBar(
    SEP, rep(list(range(SEP, na.rm=TRUE)), ncol(SEP)),
    cols=cols, drawBorders=drawBorders, main=main, yaxt=plotSampleNames,
    plotModuleNames=plotModuleNames, yaxt.line=saxt.line, maxt.line=maxt.line,
    xlab="Module summary", border.width=border.width
  )
}

#' \code{plotDataLegend}: Plot a legend for the data matrix heatmap
#' for one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotDataLegend <- function(
  data, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, palette=data.palette(), border.width=2, 
  horizontal=TRUE, legend.main="data", legend.tick.size=0.03, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input
  #-----------------------------------------------------------------------------
  if (is.null(data))
    stop("Cannot plot legend without corresponding 'data'")
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- file.path(tempdir(), paste0(".temp-objects", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  data <- unifyDS(dynamicMatLoad(data))
  correlation <- unifyDS(dynamicMatLoad(correlation))
  network <- unifyDS(dynamicMatLoad(network))
  
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
    moduleAssignments, discovery, length(correlation), names(correlation),
    ncol(correlation[[discovery]]), colnames(correlation[[discovery]])
  )
  
  if (is.null(data[[test]]))
    stop("Cannot plot module summary vector without the corresponding 'data'")
  
  #-----------------------------------------------------------------------------
  # Get the range of data matrix for the modules in the test dataset 
  #-----------------------------------------------------------------------------
  modNodes <- getModuleVarsUnsorted(moduleAssignments, modules, discovery)
  modNodes <- modNodes %sub_in% colnames(data[[test]])
  if (length(modNodes) == 0)
    stop("None of the variables composing the module are present in the test dataset")
  rg <- range(data[[test]][,modNodes])
  
  #-----------------------------------------------------------------------------
  # Plot the legend
  #-----------------------------------------------------------------------------
  emptyPlot(c(0,1), c(0,1), bty="n")
  if (all(rg < 0)) {
    addGradientLegend(
      head(palette, length(palette)/2), rg, rg, horizontal, legend.main, 
      xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size, axis.line=laxt.line,
      border.width=border.width
    )
  } else if (all(rg > 0)) {
    addGradientLegend(
      tail(palette, length(palette)/2), rg, rg, horizontal, legend.main, 
      xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size, axis.line=laxt.line,
      border.width=border.width
    )
  } else {
    plim <- c(-max(abs(rg)), max(abs(rg)))
    addGradientLegend(
      palette, plim, rg, horizontal, legend.main, xlim=c(0,1), ylim=c(0,1), 
      tick.size=legend.tick.size, axis.line=laxt.line, border.width=border.width
    )
  }
}

#' \code{plotCorrelationLegend}: Plot a legend for the correlation structure 
#' heatmap.
#' 
#' @rdname plotTopology
#' @export
plotCorrelationLegend <- function(
  palette=correlation.palette(), border.width=2, horizontal=TRUE, 
  legend.main="correlation", legend.tick.size=0.03, laxt.line=2.5, 
  cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Render legend
  #-----------------------------------------------------------------------------
  emptyPlot(c(0,1), c(0,1), bty="n")
  addGradientLegend(
    palette, c(-1,1), c(-1,1), horizontal, legend.main,
    xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size,
    axis.line=laxt.line, border.width=border.width
  )
}

#' \code{plotNetworkLegend}: Plot a legend for the network edge weights heatmap.
#' 
#' @rdname plotTopology
#' @export
plotNetworkLegend <- function(
  palette=network.palette(), border.width=2, horizontal=TRUE, 
  legend.main="network", legend.tick.size=0.03, laxt.line=2.5, 
  cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Render legend
  #-----------------------------------------------------------------------------
  emptyPlot(c(0,1), c(0,1), bty="n")
  addGradientLegend(
    palette, c(0,1), c(0,1), horizontal, legend.main,
    xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size,
    axis.line=laxt.line, border.width=border.width
  )
}
