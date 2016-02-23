#' Plot a topological feature of network module
#' 
#' Functions for plotting the topology of a network module.
#' 
#' @inheritParams common_params
#' @inheritParams par_param
#' @inheritParams orderModules_param
#' @inheritParams plot_params
#'
#' @param symmetric logical; controls whether the correlation and network 
#'  heatmaps are drawn as symmetric (square) heatmaps or asymettric triangle 
#'  heatmaps. If symmetric, then the node and module names will also be rendered
#'  on the left axis.
#' @param palette a vector of colors to use for each plot (see details).
#' @param drawBorders logical; if \code{TRUE}, borders are drawn around the bars
#'  in \code{plotContribution}, \code{plotDegree}, and
#'  \code{plotSummary}.
#' @param plotLegend logical; controls whether a legend is drawn when using
#'  \code{plotCorrelation}, \code{plotNetwork}, or \code{plotData}.
#' @param legend.main title for the legend.
#' @param main title for each plot.
#' @param legend.position the distance from the plot to start the legend, as a
#'  proportion of the plot width.
#' @param horizontal logical; controls whether the legend is rendered 
#'  horizontally or vertically when using \code{plotDataLegend},
#'  \code{plotCorrelationLegend} or \code{plotNetworkLegend}.
#'
#' @details
#'  \subsection{Input data structure:}{
#'   The \link[=modulePreservation]{preservation of network modules} in a second
#'   dataset is quantified by measuring the preservation of topological
#'   properties between the \emph{discovery} and \emph{test} datasets. These 
#'   properties are calculated not only from the interaction networks inferred
#'   in each dataset, but also from the data used to infer those networks (e.g.
#'   gene expression data) as well as the correlation structure between 
#'   variables/nodes. Thus, all functions in the \code{NetRep} package have the 
#'   following arguments: 
#'   \itemize{
#'     \item{\code{network}:}{
#'       a list of interaction networks, one for each dataset.
#'     }
#'     \item{\code{data}:}{
#'       a list of data matrices used to infer those networks, one for each 
#'       dataset.
#'     }
#'     \item{\code{correlation}:}{
#'      a list of matrices containing the pairwise correlation coefficients 
#'      between variables/nodes in each dataset.
#'     } 
#'     \item{\code{moduleAssignments}:}{
#'      a list of vectors, one for each \emph{discovery} dataset, containing 
#'      the module assignments for each node in that dataset.
#'     }
#'     \item{\code{modules}:}{
#'      a list of vectors, one for each \emph{discovery} dataset, containing
#'      the names of the modules from that dataset to analyse.  
#'     }
#'     \item{\code{discovery}:}{
#'       a vector indicating the names or indices of the previous arguments' 
#'       lists to use as the \emph{discovery} dataset(s) for the analyses.
#'     }
#'     \item{\code{test}:}{
#'       a list of vectors, one vector for each \emph{discovery} dataset, 
#'       containing the names or indices of the \code{network}, \code{data}, and 
#'       \code{correlation} argument lists to use as the \emph{test} dataset(s) 
#'       for the analysis of each \emph{discovery} dataset.
#'     }
#'   }
#'   
#'   The formatting of these arguments is not strict: each function will attempt
#'   to make sense of the user input. For example, if there is only one 
#'   \code{discovery} dataset, then input to the \code{moduleAssigments} and 
#'   \code{test} arguments may be vectors, rather than lists. If the node and
#'   sample ordering is being calculated within the same dataset being 
#'   visualised, then the \code{discovery} and \code{test} arguments do
#'   not need to be specified, and the input matrices for the \code{network},
#'   \code{data}, and \code{correlation} arguments do not need to be wrapped in
#'   a list.
#' }
#' \subsection{'bigMatrix' input data:}{
#'   Although the \code{data}, \code{correlation}, and \code{network} arguments 
#'   expect data to be provided in the \code{\link{bigMatrix}}, they can be 
#'   provided as regular \code{\link[base]{matrix}} objects, in which case they 
#'   will be temporarily converted to \code{bigMatrix} objects. \strong{This is 
#'   not recommended}, as each matrix will be copied into memory of each 
#'   parallel R session, resulting in much higher memory usage (\code{bigMatrix}
#'   objects can be simultaneously accessed from multiple R sessions with no 
#'   additional memory overhead), and increased computation time for the 
#'   conversion and copying processes. It is therefore strongly recommended that
#'   the user save their data separately as \code{\link{bigMatrix}} objects
#'   prior to running the permutation procedure or using any other package
#'   function. This is also useful for other analyses, as \code{bigMatrix} 
#'   objects can be instantaneously loaded into any future R session.
#'   
#'   Alternatively, the \code{data}, \code{correlation}, and \code{network} 
#'   arguments will also accept file paths to tabular data or \code{bigMatrix}
#'   backingfiles.
#' }
#' \subsection{Memory usage:}{
#'   Provided there are no additional objects in the R session the
#'   permutation procedure will use only the memory required to store each 
#'   matrix once, along with an additional 200 MB per core used by each vanilla
#'   R session.
#' }
#' \subsection{Node, sample, and module ordering:}{
#'   By default, nodes are ordered in decreasing order of \emph{weighted degree}
#'   in the \code{discovery} dataset (see \code{\link{nodeOrder}}). 
#'   This facilitates the visual comparison of modules across datasets, as the 
#'   node ordering will be preserved. Missing nodes are colored in grey. If
#'   \code{orderNodesBy} is "test" nodes will instead be ordered by 
#'   weighted degree in the \code{test} dataset. If "none" nodes are 
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
#' \subsection{Normalised degree:}{
#'   The weighted degree is normalised by the maximum connectivity in
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
#'   the bars in \code{plotDegree}, \code{plotContribution}, and 
#'   \code{plotSummary}.
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
#'   \code{plotDegree} expects \code{palette} to be a single color, a 
#'   vector of colors, one for each node, or a vector of colors to be repeated.
#'   
#'   \code{plotContribution} expects \code{palette} to be a vector 
#'   containing two colors, the first to be used for nodes with negative node
#'   contribution values, and the second to be used for nodes with positive node
#'   contribution values. 
#'   
#'   \code{plotData} and \code{plotDataLegend} expect the \code{palette}
#'   argument to be a vector of colors to interpolate over when plotting the
#'   'data.' In order to accomodate data matrices with different ranges,
#'   these functions expect the palette to be a diverging set of colors with a
#'   centre value of 0 (e.g. white). The colors drawn are balanced around 0: 
#'   i.e. positive and negative values of data will have the same intensity on a 
#'   diverging color palette regardless of the actual range of the data.
#'   
#'   \code{plotSummary} expects \code{palette} to be a vector 
#'   containing two colors, the first to be used for nodes with a negative 
#'   module summary measurement, and the second to be used for genes with a 
#'   positive module summary measurment, regardless of whether the data matrix
#'   is centred around 0.
#' }
#' 
#' @seealso
#' \code{\link{plotModule}} for a combined plot showing all topological 
#' properties for a network module.
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
#' # Plot the data for module 1, 2 and 4 in the discovery dataset
#' plotData(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4)
#' )
#' 
#' # Symmetric = TRUE gives a traditional heatmap for the correlation structure
#' # and weighted network
#' plotCorrelation(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4), symmetric=TRUE
#' )
#' 
#' # While the default is to render only one half of the (symmetric) matrix
#' plotNetwork(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4)
#' )
#' 
#' # Plot the degree of nodes in each module in the test dataset, but show them
#' # in the same order as the discovery dataset to compare how node degree 
#' # changes
#' plotDegree(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4), discovery="discovery",
#'   test="test"
#' )
#' 
#' # Alternatively nodes can be order on the plot by degree in the test dataset
#' plotDegree(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4), discovery="discovery",
#'   test="test", orderNodesBy="test"
#' )
#' 
#' # Arbitrary subsets can be plotted:
#' plotContribution(
#'   data=data_list[[1]][, 1:10], correlation=correlation_list[[1]][1:10, 1:10], 
#'   network=network_list[[1]][1:10, 1:10], orderNodesBy="none" 
#' )
#' 
#' # Plot the module summary vectors for multiple modules:
#' plotSummary(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4), discovery="discovery",
#'   test="test", orderSamplesBy="test"
#' )
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
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, verbose=TRUE,
  orderSamplesBy="test", orderNodesBy="discovery",
  orderModules=TRUE, plotNodeNames=TRUE, plotSampleNames=TRUE, plotModuleNames,
  main="", palette=data.palette(), border.width=2, plotLegend=TRUE, 
  legend.main="Data", gaxt.line=-0.5, saxt.line=-0.5, maxt.line=3, 
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  if (is.null(data))
    stop("Cannot plot data matrix without 'data'")
  
  # Check plot-specific arguments
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
  if (!is.logical(orderModules) || is.na(orderModules) || length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # At this time, we can only plot within one dataset.
  if ((!is.null(discovery) && (!is.vector(discovery) || length(discovery) > 1)) ||
      (!is.null(test) && (!is.vector(test) || length(test) > 1))) {
    stop("only 1 'discovery' and 'test' dataset can be specified when plotting")
  }
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # set up 'discovery' as 'test' so we can use it on 'netPropsInternal'
  discAsTest <- list(discovery)
  names(discAsTest) <- discovery
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (missing(plotModuleNames))
    plotModuleNames <- length(mods) > 1
  
  if (is.null(data[[ti]]))
    stop("Cannot plot data matrix without 'data'")
  
  if ((orderSamplesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderSamplesBy' dataset") 
  }
  
  if (orderSamplesBy == "discovery" && 
      sum(rownames(scaledData[di]) %in% rownames(scaledData[[ti]])) == 0) {
    stop("'orderBySamples' can only be ", '"discovery"', " when the same",
         " samples are present in both the 'discovery' and 'test' datasets")
  }
  
  if ((orderModules && length(mods) > 1) && 
      (orderNodesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderNodesBy' dataset and ",
         "'orderModules' = 'TRUE'") 
  }
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and samples in the 'test' dataset by the dataset 
  # specified in 'orderNodesBy' and 'orderSamplesBy'.
  #-----------------------------------------------------------------------------
  # Calculate the network properties in the dataset we're plotting.
  testProps <- netPropsInternal(
    scaledData, correlation, network, moduleAssignments, 
    modules, discovery, test, datasetNames, FALSE
  )
  
  # Case 1: we want to order nodes by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderNodesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    discProps <- netPropsInternal(
      scaledData, correlation, network, moduleAssignments, 
      modules, discovery, discAsTest, datasetNames, FALSE
    )
    nodeOrder <- nodeOrderInternal(
      discProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE
    )[[di]][[di]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  } 
  # Case 2: order nodes as they're provided by the user
  else if (orderNodesBy == "none") {
    moduleOrder <- names(testProps[[di]][[ti]])
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(testProps[[di]][[ti]][[mi]]$degree)
    }
  } 
  # Case 3: order nodes by their degree in the test network.
  else {
    # Order modules and samples by the test network
    nodeOrder <- nodeOrderInternal(
      testProps, orderModules, simplify=FALSE, verbose, na.rm=TRUE
    )[[di]][[ti]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  }
  
  # Case 1: we want to order samples by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderSamplesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    if (!exists("discProps")) {
      discProps <- netPropsInternal(
        scaledData, correlation, network, moduleAssignments, 
        modules, discovery, discAsTest, datasetNames, FALSE
      )
    }
    sampleOrder <- sampleOrderInternal(discProps, verbose, FALSE)
    sampleOrder <- sampleOrder[[di]][[di]][[moduleOrder[1]]]
  } 
  # Case 2: order samples as they're provided by the user
  else if (orderSamplesBy == "none") {
    sampleOrder <- seq_along(testProps[[di]][[ti]][[moduleOrder[1]]]$summary)
  } 
  # Case 3: order samples by their degree in the test network.
  else {
    # Order modules and samples by the test network
    sampleOrder <- sampleOrderInternal(testProps, verbose, TRUE)
    sampleOrder <- sampleOrder[[di]][[ti]][[moduleOrder[1]]]
  }
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  
  # Case 1: di == ti. Plotting within the same dataset => nothing missing.
  # Case 2: orderBy == di: those missing in the discovery should have grey bars.
  # Case 3: orderBy == ti: ordering within the same dataset => nothing missing.
  
  if (orderNodesBy == "discovery" && di != ti) {
    na.pos.x <- which(nodeOrder %nin% colnames(network[[ti]]))
    if (length(na.pos.x) > 0) {
      presentNodes <- nodeOrder[-na.pos.x]
    } else {
      presentNodes <- nodeOrder
    }
  } else {
    na.pos.x <- vector()
    presentNodes <- nodeOrder
  }
  
  if (orderSamplesBy == "discovery" && di != ti) {
    na.pos.y <- which(sampleOrder %nin% colnames(network[[ti]]))
    if (length(na.pos.y) > 0) {
      presentSamples <- sampleOrder[-na.pos.y]
    } else {
      presentSamples <- sampleOrder
    }
  } else {
    na.pos.y <- vector()
    presentSamples <- sampleOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the data matrix
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  # First we need to set up the color palette for the data, which 
  # includes the fact that the range may be unbalanced around 0.
  dat <- data[[ti]][presentSamples, presentNodes]
  range.dat <- range(dat)
  if (all(range.dat > 0)) {
    palette <- tail(data.palette(), length(data.palette())/2)
    range.pal <- range.dat
  } else if (all(range.dat < 0)) {
    palette <- head(data.palette(), length(data.palette())/2)
    range.pal <- range.dat
  } else {
    palette <- data.palette()
    range.pal <- c(-max(abs(range.dat)), max(abs(range.dat)))
  }
  
  # Axis tick labels
  xaxt <- NULL
  if (plotNodeNames)
    xaxt <- nodeOrder
  
  yaxt <- NULL
  if (plotSampleNames)
    yaxt <- sampleOrder
  
  # Add space for the legend
  if (plotLegend) {
    old.mar <- par("mar")
    par(mar=par("mar")+c(0,0,0,2))
    on.exit({
      par(mar=old.mar)
    }, add=TRUE)
  }
    
  # Plot
  plotSquareHeatmap(
    dat, palette, vlim=range.pal, legend.lim=range.dat,
    moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.y, 
    xaxt=xaxt, yaxt=yaxt, plotLegend=plotLegend, main=main,
    legend.main=legend.main, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, yaxt.line=saxt.line, legend.tick.size=legend.tick.size,
    laxt.line=laxt.line, legend.line=legend.position, maxt.line=maxt.line,
    border.width=border.width
  )
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotCorrelation}: Plot a heatmap of the correlation structure for one 
#' or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotCorrelation <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, verbose=TRUE,
  orderNodesBy="discovery", symmetric=FALSE, orderModules=TRUE, 
  plotNodeNames=TRUE, plotModuleNames, main="", palette=correlation.palette(), 
  border.width=2, plotLegend=TRUE, legend.main="Correlation", gaxt.line=-0.5, 
  maxt.line=3, legend.position, legend.tick.size, laxt.line, cex.axis=0.8, 
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  # Check plot-specific arguments
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]

  if (!is.logical(orderModules) || is.na(orderModules) || length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # At this time, we can only plot within one dataset.
  if ((!is.null(discovery) && (!is.vector(discovery) || length(discovery) > 1)) ||
      (!is.null(test) && (!is.vector(test) || length(test) > 1))) {
    stop("only 1 'discovery' and 'test' dataset can be specified when plotting")
  }
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # set up 'discovery' as 'test' so we can use it on 'netPropsInternal'
  discAsTest <- list(discovery)
  names(discAsTest) <- discovery
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (missing(plotModuleNames))
    plotModuleNames <- length(mods) > 1
  
  if ((orderModules && length(mods) > 1) && 
      (orderNodesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderNodesBy' dataset and ",
         "'orderModules' = 'TRUE'") 
  }
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes in the 'test' dataset by the dataset specified in 
  # 'orderNodesBy'.
  #-----------------------------------------------------------------------------

  # Calculate the network properties in the dataset we're plotting.
  testProps <- netPropsInternal(
    scaledData, correlation, network, moduleAssignments, 
    modules, discovery, test, datasetNames, FALSE
  )
  
  # Case 1: we want to order nodes by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderNodesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    discProps <- netPropsInternal(
      scaledData, correlation, network, moduleAssignments, 
      modules, discovery, discAsTest, datasetNames, FALSE
    )
    nodeOrder <- nodeOrderInternal(
      discProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE
    )[[di]][[di]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  } 
  # Case 2: order nodes as they're provided by the user
  else if (orderNodesBy == "none") {
    moduleOrder <- names(testProps[[di]][[ti]])
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(testProps[[di]][[ti]][[mi]]$degree)
    }
  } 
  # Case 3: order nodes by their degree in the test network.
  else {
    # Order modules and samples by the test network
    nodeOrder <- nodeOrderInternal(
      testProps, orderModules, simplify=FALSE, verbose, na.rm=TRUE
    )[[di]][[ti]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes from the 'discovery' dataset not present in the 'test' 
  # dataset.
  #-----------------------------------------------------------------------------
  
  # Case 1: di == ti. Plotting within the same dataset => nothing missing.
  # Case 2: orderBy == di: those missing in the discovery should have grey bars.
  # Case 3: orderBy == ti: ordering within the same dataset => nothing missing.
  
  if (orderNodesBy == "discovery" && di != ti) {
    na.pos.x <- which(nodeOrder %nin% colnames(network[[ti]]))
    if (length(na.pos.x) > 0) {
      presentNodes <- nodeOrder[-na.pos.x]
    } else {
      presentNodes <- nodeOrder
    }
  } else {
    na.pos.x <- vector()
    presentNodes <- nodeOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the correlation structure
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  # Plot axis tick labels?
  gaxt <- NULL
  if (plotNodeNames)
    gaxt <- nodeOrder
  
  if (symmetric) {
    # Set defaults
    if (missing(legend.position))
      legend.position <- 0.15
    if (missing(laxt.line))
      laxt.line <- 3
    if (missing(legend.tick.size))
      legend.tick.size <- 0.03
    
    # Add space for the legend
    if (plotLegend) {
      old.mar <- par("mar")
      par(mar=par("mar")+c(0,0,0,3))
      on.exit({
        par(mar=old.mar)
      }, add=TRUE)
    }
    
    plotSquareHeatmap(
      correlation[[ti]][presentNodes, presentNodes], palette, c(-1, 1), 
      moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.x, 
      xaxt=gaxt, yaxt=gaxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, yaxt.line=gaxt.line, border.width=border.width,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line
    )
  } else {
    # Set defaults 
    if (missing(legend.position))
      legend.position <- 0.1
    if (missing(laxt.line))
      laxt.line <- 2.5
    if (missing(legend.tick.size))
      legend.tick.size <- 0.025
    
    plotTriangleHeatmap(
      correlation[[ti]][presentNodes, presentNodes], palette, c(-1, 1),
      moduleAssignments[[di]][nodeOrder], na.pos.x, xaxt=gaxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=gaxt.line,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line, 
      border.width=border.width
    )
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotNetwork}: Plot a heatmap of the edge weights for one or more
#' network modules.
#' 
#' @rdname plotTopology
#' @export
plotNetwork <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, verbose=TRUE,
  orderNodesBy="discovery", symmetric=FALSE, orderModules=TRUE, 
  plotNodeNames=TRUE, plotModuleNames, main="", palette=network.palette(), 
  border.width=2, plotLegend=TRUE, legend.main="Edge weight", gaxt.line=-0.5, 
  maxt.line=3, legend.position, legend.tick.size, laxt.line, cex.axis=0.8,
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  # Check plot-specific arguments
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  
  if (!is.logical(orderModules) || is.na(orderModules) || length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # At this time, we can only plot within one dataset.
  if ((!is.null(discovery) && (!is.vector(discovery) || length(discovery) > 1)) ||
      (!is.null(test) && (!is.vector(test) || length(test) > 1))) {
    stop("only 1 'discovery' and 'test' dataset can be specified when plotting")
  }
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # set up 'discovery' as 'test' so we can use it on 'netPropsInternal'
  discAsTest <- list(discovery)
  names(discAsTest) <- discovery
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (missing(plotModuleNames))
    plotModuleNames <- length(mods) > 1
  
  if ((orderModules && length(mods) > 1) && 
      (orderNodesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderNodesBy' dataset and ",
         "'orderModules' = 'TRUE'") 
  }
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes in the 'test' dataset by the dataset specified in 
  # 'orderNodesBy'.
  #-----------------------------------------------------------------------------
  
  # Calculate the network properties in the dataset we're plotting.
  testProps <- netPropsInternal(
    scaledData, correlation, network, moduleAssignments, 
    modules, discovery, test, datasetNames, FALSE
  )
  
  # Case 1: we want to order nodes by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderNodesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    discProps <- netPropsInternal(
      scaledData, correlation, network, moduleAssignments, 
      modules, discovery, discAsTest, datasetNames, FALSE
    )
    nodeOrder <- nodeOrderInternal(
      discProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE
    )[[di]][[di]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  } 
  # Case 2: order nodes as they're provided by the user
  else if (orderNodesBy == "none") {
    moduleOrder <- names(testProps[[di]][[ti]])
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(testProps[[di]][[ti]][[mi]]$degree)
    }
  } 
  # Case 3: order nodes by their degree in the test network.
  else {
    # Order modules and samples by the test network
    nodeOrder <- nodeOrderInternal(
      testProps, orderModules, simplify=FALSE, verbose, na.rm=TRUE
    )[[di]][[ti]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes from the 'discovery' dataset not present in the 'test' 
  # dataset.
  #-----------------------------------------------------------------------------
  
  # Case 1: di == ti. Plotting within the same dataset => nothing missing.
  # Case 2: orderBy == di: those missing in the discovery should have grey bars.
  # Case 3: orderBy == ti: ordering within the same dataset => nothing missing.
  
  if (orderNodesBy == "discovery" && di != ti) {
    na.pos.x <- which(nodeOrder %nin% colnames(network[[ti]]))
    if (length(na.pos.x) > 0) {
      presentNodes <- nodeOrder[-na.pos.x]
    } else {
      presentNodes <- nodeOrder
    }
  } else {
    na.pos.x <- vector()
    presentNodes <- nodeOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the network edge weights
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  # Plot axis tick labels?
  gaxt <- NULL
  if (plotNodeNames)
    gaxt <- nodeOrder
  
  if (symmetric) {
    # Set defaults
    if (missing(legend.position))
      legend.position <- 0.15
    if (missing(laxt.line))
      laxt.line <- 3
    if (missing(legend.tick.size))
      legend.tick.size <- 0.03
    
    # Add space for the legend
    if (plotLegend) {
      old.mar <- par("mar")
      par(mar=par("mar")+c(0,0,0,3))
      on.exit({
        par(mar=old.mar)
      }, add=TRUE)
    }
    
    plotSquareHeatmap(
      network[[ti]][presentNodes, presentNodes], palette, c(0, 1), 
      moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.x, 
      xaxt=gaxt, yaxt=gaxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, yaxt.line=gaxt.line, 
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line,
      border.width=border.width
    )
  } else {
    # Set defaults 
    if (missing(legend.position))
      legend.position <- 0.1
    if (missing(laxt.line))
      laxt.line <- 2.5
    if (missing(legend.tick.size))
      legend.tick.size <- 0.025
    
    plotTriangleHeatmap(
      network[[ti]][presentNodes, presentNodes], palette, c(0, 1), 
      moduleAssignments[[di]][nodeOrder], na.pos.x,xaxt=gaxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=gaxt.line,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line, 
      border.width=border.width
    )
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotContribution}: Plot a bar chart of the module membership for
#' one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotContribution <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, verbose=TRUE,
  orderNodesBy="discovery", orderModules=TRUE, plotNodeNames=TRUE, 
  plotModuleNames, main="", border.width=2, palette=c("#313695", "#a50026"), 
  drawBorders=FALSE, gaxt.line=-0.5, maxt.line=3, cex.axis=0.8, cex.lab=1, 
  cex.main=1.2
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  if (is.null(data))
    stop("Cannot plot node contribution without 'data'")
  
  # Check plot-specific arguments
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]

  if (!is.logical(orderModules) || is.na(orderModules) || length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # At this time, we can only plot within one dataset.
  if ((!is.null(discovery) && (!is.vector(discovery) || length(discovery) > 1)) ||
      (!is.null(test) && (!is.vector(test) || length(test) > 1))) {
    stop("only 1 'discovery' and 'test' dataset can be specified when plotting")
  }
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # set up 'discovery' as 'test' so we can use it on 'netPropsInternal'
  discAsTest <- list(discovery)
  names(discAsTest) <- discovery
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (missing(plotModuleNames))
    plotModuleNames <- length(mods) > 1
  
  if (is.null(data[[ti]]))
    stop("Cannot plot node contribution without 'data'")
  
  if ((orderModules && length(mods) > 1) && 
      (orderNodesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderNodesBy' dataset and ",
         "'orderModules' = 'TRUE'") 
  }
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and samples in the 'test' dataset by the dataset 
  # specified in 'orderNodesBy'.
  #-----------------------------------------------------------------------------
  # Calculate the network properties in the dataset we're plotting.
  testProps <- netPropsInternal(
    scaledData, correlation, network, moduleAssignments, 
    modules, discovery, test, datasetNames, FALSE
  )
  
  # Case 1: we want to order nodes by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderNodesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    discProps <- netPropsInternal(
      scaledData, correlation, network, moduleAssignments, 
      modules, discovery, discAsTest, datasetNames, FALSE
    )
    nodeOrder <- nodeOrderInternal(
      discProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE
    )[[di]][[di]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  } 
  # Case 2: order nodes as they're provided by the user
  else if (orderNodesBy == "none") {
    moduleOrder <- names(testProps[[di]][[ti]])
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(testProps[[di]][[ti]][[mi]]$degree)
    }
  } 
  # Case 3: order nodes by their degree in the test network.
  else {
    # Order modules and samples by the test network
    nodeOrder <- nodeOrderInternal(
      testProps, orderModules, simplify=FALSE, verbose, na.rm=TRUE
    )[[di]][[ti]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Node contribution
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  testProps <- testProps[[di]][[ti]]
  
  # node contribution
  nodeContribVec <- foreach(mi = seq_along(testProps), .combine=c) %do% {
    testProps[[mi]]$contribution
  }
  nodeContribVec <- nodeContribVec[nodeOrder]
  
  # Plot bar chart
  plotBar(
    nodeContribVec, c(-1,1), moduleAssignments[[di]][nodeOrder],
    ifelse(nodeContribVec > 0, palette[2], palette[1]), drawBorders=drawBorders,
    xaxt=plotNodeNames, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, maxt.line=maxt.line, main=main,
    ylab="Node contribution", border.width=border.width
  )
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotDegree:} Plot a bar chart of the normalised intramodular 
#' connectivity (see details) for one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotDegree <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, verbose=TRUE,
  orderNodesBy="discovery", orderModules=TRUE, plotNodeNames=TRUE, 
  plotModuleNames, main="", palette="#feb24c", border.width=2, 
  drawBorders=FALSE, gaxt.line=-0.5, maxt.line=3, cex.axis=0.8, cex.lab=1, 
  cex.main=1.2
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  # Check plot-specific arguments
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]

  if (!is.logical(orderModules) || is.na(orderModules) || length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # At this time, we can only plot within one dataset.
  if ((!is.null(discovery) && (!is.vector(discovery) || length(discovery) > 1)) ||
      (!is.null(test) && (!is.vector(test) || length(test) > 1))) {
    stop("only 1 'discovery' and 'test' dataset can be specified when plotting")
  }
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # set up 'discovery' as 'test' so we can use it on 'netPropsInternal'
  discAsTest <- list(discovery)
  names(discAsTest) <- discovery
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (missing(plotModuleNames))
    plotModuleNames <- length(mods) > 1

  if ((orderModules && length(mods) > 1) && 
      (orderNodesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderNodesBy' dataset and ",
         "'orderModules' = 'TRUE'") 
  }
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and samples in the 'test' dataset by the dataset 
  # specified in 'orderNodesBy'.
  #-----------------------------------------------------------------------------
  # Calculate the network properties in the dataset we're plotting.
  testProps <- netPropsInternal(
    scaledData, correlation, network, moduleAssignments, 
    modules, discovery, test, datasetNames, FALSE
  )
  
  # Case 1: we want to order nodes by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderNodesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    discProps <- netPropsInternal(
      scaledData, correlation, network, moduleAssignments, 
      modules, discovery, discAsTest, datasetNames, FALSE
    )
    nodeOrder <- nodeOrderInternal(
      discProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE
    )[[di]][[di]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  } 
  # Case 2: order nodes as they're provided by the user
  else if (orderNodesBy == "none") {
    moduleOrder <- names(testProps[[di]][[ti]])
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(testProps[[di]][[ti]][[mi]]$degree)
    }
  } 
  # Case 3: order nodes by their degree in the test network.
  else {
    # Order modules and samples by the test network
    nodeOrder <- nodeOrderInternal(
      testProps, orderModules, simplify=FALSE, verbose, na.rm=TRUE
    )[[di]][[ti]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Node contribution
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  testProps <- testProps[[di]][[ti]]
  
  # (Normalised) weighted degree vector
  wDegreeVec <- foreach(mi = seq_along(testProps), .combine=c) %do% {
    testProps[[mi]]$degree/max(na.omit(testProps[[mi]]$degree))
  }
  wDegreeVec <- wDegreeVec[nodeOrder]
  
  # Plot bar chart
  plotBar(
    wDegreeVec, c(0,1), moduleAssignments[[di]][nodeOrder],
    palette, drawBorders=drawBorders,
    xaxt=plotNodeNames, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, maxt.line=maxt.line, main=main,
    ylab="Weighted degree", border.width=border.width
  )
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotSummary}: Plot bar charts of the module summary vectors of 
#' one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotSummary <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, verbose=TRUE,
  orderSamplesBy="test", orderNodesBy="discovery", orderModules=TRUE, 
  plotSampleNames=TRUE, plotModuleNames, main="", 
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  if (is.null(data))
    stop("Cannot plot module summary without 'data'")
  
  # Check plot-specific arguments
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderNodesBy <- orderByArgs[pmatch(orderNodesBy, orderByArgs, nomatch=3)]
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
  if (!is.logical(orderModules) || is.na(orderModules) || length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # At this time, we can only plot within one dataset.
  if ((!is.null(discovery) && (!is.vector(discovery) || length(discovery) > 1)) ||
      (!is.null(test) && (!is.vector(test) || length(test) > 1))) {
    stop("only 1 'discovery' and 'test' dataset can be specified when plotting")
  }
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # set up 'discovery' as 'test' so we can use it on 'netPropsInternal'
  discAsTest <- list(discovery)
  names(discAsTest) <- discovery
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (missing(plotModuleNames))
    plotModuleNames <- length(mods) > 1
  
  if (is.null(data[[ti]]))
    stop("Cannot plot module summary without 'data'")
  
  if ((orderSamplesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderSamplesBy' dataset") 
  }
  
  if (orderSamplesBy == "discovery" && 
      sum(rownames(scaledData[di]) %in% rownames(scaledData[[ti]])) == 0) {
    stop("'orderBySamples' can only be ", '"discovery"', " when the same",
         " samples are present in both the 'discovery' and 'test' datasets")
  }
  
  if ((orderModules && length(mods) > 1) && 
      (orderNodesBy == "discovery" && is.null(scaledData[[di]]))) {
    stop("'data' not provided for 'orderNodesBy' dataset and ",
         "'orderModules' = 'TRUE'") 
  }
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and samples in the 'test' dataset by the dataset 
  # specified in 'orderNodesBy' and 'orderSamplesBy'.
  #-----------------------------------------------------------------------------
  # Calculate the network properties in the dataset we're plotting.
  testProps <- netPropsInternal(
    scaledData, correlation, network, moduleAssignments, 
    modules, discovery, test, datasetNames, FALSE
  )
  
  # Case 1: we want to order nodes by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderNodesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    discProps <- netPropsInternal(
      scaledData, correlation, network, moduleAssignments, 
      modules, discovery, discAsTest, datasetNames, FALSE
    )
    nodeOrder <- nodeOrderInternal(
      discProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE
    )[[di]][[di]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  } 
  # Case 2: order nodes as they're provided by the user
  else if (orderNodesBy == "none") {
    moduleOrder <- names(testProps[[di]][[ti]])
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(testProps[[di]][[ti]][[mi]]$degree)
    }
  } 
  # Case 3: order nodes by their degree in the test network.
  else {
    # Order modules and samples by the test network
    nodeOrder <- nodeOrderInternal(
      testProps, orderModules, simplify=FALSE, verbose, na.rm=TRUE
    )[[di]][[ti]]
    moduleOrder <- names(nodeOrder)
    nodeOrder <- unlist(nodeOrder)
  }
  
  # Case 1: we want to order samples by the discovery dataset, which if different
  # to the test dataset, we need to recalculate the weighted degree for the 
  # node order.
  if (orderSamplesBy == "discovery" && di != ti) {
    # This skips all of the data verification
    if (!exists("discProps")) {
      discProps <- netPropsInternal(
        scaledData, correlation, network, moduleAssignments, 
        modules, discovery, discAsTest, datasetNames, FALSE
      )
    }
    sampleOrder <- sampleOrderInternal(discProps, verbose, FALSE)
    sampleOrder <- sampleOrder[[di]][[di]][[moduleOrder[1]]]
  } 
  # Case 2: order samples as they're provided by the user
  else if (orderSamplesBy == "none") {
    sampleOrder <- seq_along(testProps[[di]][[ti]][[moduleOrder[1]]]$summary)
  } 
  # Case 3: order samples by their degree in the test network.
  else {
    # Order modules and samples by the test network
    sampleOrder <- sampleOrderInternal(testProps, verbose, TRUE)
    sampleOrder <- sampleOrder[[di]][[ti]][[moduleOrder[1]]]
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  
  # Case 1: di == ti. Plotting within the same dataset => nothing missing.
  # Case 2: orderBy == di: those missing in the discovery should have grey bars.
  # Case 3: orderBy == ti: ordering within the same dataset => nothing missing.
  
  if (orderSamplesBy == "discovery" && di != ti) {
    na.pos.y <- which(sampleOrder %nin% colnames(network[[ti]]))
    if (length(na.pos.y) > 0) {
      presentSamples <- sampleOrder[-na.pos.y]
    } else {
      presentSamples <- sampleOrder
    }
  } else {
    na.pos.y <- vector()
    presentSamples <- sampleOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Node contribution
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  testProps <- testProps[[di]][[ti]]
  
  # Summary profile matrix
  summaries <- foreach(mi = moduleOrder, .combine=cbind) %do% {
    matrix(
      insert.nas(testProps[[mi]]$summary[presentSamples], na.pos.y),
      ncol=1
    )
  }
  colnames(summaries) <- moduleOrder
  rownames(summaries) <- sampleOrder
  
  # and respective colors
  summaries.cols <- matrix(data.palette()[1], nrow(summaries), ncol(summaries))
  summaries.cols[summaries > 0] <- tail(data.palette(), 1)
  
  # Plot bar chart
  plotMultiBar(
    summaries, rep(list(range(summaries, na.rm=TRUE)), ncol(summaries)),
    cols=summaries.cols, drawBorders=drawBorders, main=main, 
    yaxt=plotSampleNames, plotModuleNames=plotModuleNames, yaxt.line=saxt.line, 
    maxt.line=maxt.line, xlab="Module summary", border.width=border.width
  )
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotDataLegend}: Plot a legend for the data matrix heatmap
#' for one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotDataLegend <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, verbose=TRUE, 
  palette=data.palette(), border.width=2, horizontal=TRUE, legend.main="Data", 
  legend.tick.size=0.03, laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  if (is.null(data))
    stop("Cannot plot data legend without 'data'")
  
  # Check plot-specific arguments
  if (class(legend.main) != "character")
    stop("'legend.main' must be a characer vector")
  
  # At this time, we can only plot within one dataset.
  if ((!is.null(discovery) && (!is.vector(discovery) || length(discovery) > 1)) ||
      (!is.null(test) && (!is.vector(test) || length(test) > 1))) {
    stop("only 1 'discovery' and 'test' dataset can be specified when plotting")
  }
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  scaledData <- finput$scaledData
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]

  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (is.null(data[[ti]]))
    stop("Cannot plot data legend without 'data'")
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get the range of data matrix for the modules in the test dataset 
  #-----------------------------------------------------------------------------
  modNodes <- getModuleVarsUnsorted(moduleAssignments, mods, di)
  modNodes <- modNodes %sub_in% colnames(data[[ti]])
  if (length(modNodes) == 0)
    stop("None of the variables composing the module are present in the test dataset")
  rg <- range(data[[ti]][,modNodes])
  
  #-----------------------------------------------------------------------------
  # Plot the legend
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
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
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotCorrelationLegend}: Plot a legend for the correlation structure 
#' heatmap.
#' 
#' @rdname plotTopology
#' @export
plotCorrelationLegend <- function(
  palette=correlation.palette(), border.width=2, horizontal=TRUE, 
  legend.main="correlation", legend.tick.size=0.03, laxt.line=2.5, 
  cex.axis=0.8, cex.lab=1, cex.main=1.2, verbose=TRUE
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
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotNetworkLegend}: Plot a legend for the network edge weights heatmap.
#' 
#' @rdname plotTopology
#' @export
plotNetworkLegend <- function(
  palette=network.palette(), border.width=2, horizontal=TRUE, 
  legend.main="network", legend.tick.size=0.03, laxt.line=2.5, 
  cex.axis=0.8, cex.lab=1, cex.main=1.2, verbose=TRUE
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
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}
