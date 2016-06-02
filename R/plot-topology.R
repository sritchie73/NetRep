#' Plot a topological feature of network module
#' 
#' Functions for plotting the topology of a network module.
#' 
#' @inheritParams common_params
#' @inheritParams orderModules_param
#' @inheritParams plot_params
#'
#' @param symmetric logical; controls whether the correlation and network 
#'  heatmaps are drawn as symmetric (square) heatmaps or asymettric triangle 
#'  heatmaps. If symmetric, then the node and module names will also be rendered
#'  on the left axis.
#' @param plotLegend logical; controls whether a legend is drawn when using
#'  \code{plotCorrelation}, \code{plotNetwork}, or \code{plotData}.
#' @param legend.main title for the legend.
#' @param legend.position the distance from the plot to start the legend, as a
#'  proportion of the plot width.
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
#'   parallel R session, resulting in much higher memory usage (see section on
#'   Memory usage ), and increased computation time for the conversion and
#'   copying processes. It is therefore strongly recommended that the user save
#'   their data separately as \code{\link{bigMatrix}} objects prior to running
#'   the permutation procedure or using any other package function. This is also
#'   useful for other analyses, as \code{bigMatrix} objects can be
#'   instantaneously loaded into any future R session.
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
#'   in the \code{discovery} dataset (see \code{\link{nodeOrder}}). Missing 
#'   nodes are colored in grey. This facilitates the visual comparison of 
#'   modules across datasets, as the node ordering will be preserved. 
#'   
#'   Alternatively, a vector containing the names or indices of one or more
#'   datasets can be provided to the \code{orderNodesBy} argument. If a single
#'   dataset is provided, then nodes will be ordered in decreasing order of
#'   \emph{weighted degree} in that dataset. If multiple datasets are provided
#'   then the \emph{weighted degree} will be calculated in each dataset and
#'   nodes will be ordered in decreasing order of average weighted degree. This
#'   is useful for obtaining a robust ordering of nodes by importance across
#'   independent datasets where a module is reproducible.
#'   
#'   Ordering of nodes by \emph{weighted degree} can be suppressed by setting
#'   \code{orderNodesBy} to \code{NA}, in which case nodes will be ordered as in
#'   the matrices provided in the \code{data}, \code{correlation} and
#'   \code{network} arguments.
#'   
#'   When multiple modules are drawn, modules are ordered by the similarity
#'   of their summary vectors in the dataset(s) specified in \code{orderNodesBy}
#'   argument. If multiple datasets are provided to the \code{orderNodesBy}
#'   argument then the module summary vectors are concatenated across datasets.
#'   
#'   By default, samples in the data heatmap drawn by \code{plotData} and the 
#'   bar plot drawn by \code{plotModuleSummary} are ordered in descending order
#'   of \emph{module summary} in the drawn dataset (specified by the \code{test}
#'   argument). If multiple modules are drawn, samples are ordered as per the
#'   left-most module on the plot.
#'   
#'   Alternatively, a vector containing the name or index of another dataset may
#'   be provided to the \code{orderSamplesBy} argument. In this case, samples
#'   will be ordered in descending order of \emph{module summary} in the 
#'   specified dataset. This is useful when comparing different measurements 
#'   across samples, for example, gene expression data obtained from multiple 
#'   tissues samples across the same individuals. Samples that are not present 
#'   in the drawn dataset will be colored in grey.
#'   
#'   Order of samples by \emph{module summary} can be suppressed by setting 
#'   \code{orderSamplesBy} to \code{NA}, in which case samples will be order as
#'   in the matrix provided to the \code{data} argument for the drawn dataset.
#' }
#' \subsection{Normalised degree:}{
#'   The weighted degree is normalised by the maximum weighted degree in
#'   any given module when rendered on the bar plot. This facilitates visual 
#'   comparison of multiple modules with differing sizes or densities.
#' }
#' \subsection{Plot layout and device size}{
#'   Although reasonable default values for most parameters have been provided,
#'   the rendering of axes and titles may need adjusting depending on the size
#'   of the plot window. The \code{dryRun} argument is useful for quickly 
#'   determining whether the plot will render correctly. When \code{dryRun} 
#'   is \code{TRUE} only the axes, legends, labels, and title will be drawn, 
#'   allowing for quick iteration of customisable parameters to get the plot 
#'   layout correct. 
#'   
#'   \strong{Warning}: PDF and other vectorized devices should not be used when
#'   plotting the heatmaps with more than a hundred nodes. Large files will be
#'   generated which may cause image editing programs such as Inkscape or
#'   Illustrator to crash when polishing figures for publication.
#'   
#'   If axis labels or legends are drawn off screen then the margins of the 
#'   plot should be adjusted prior to plotting using the 
#'   \code{\link[graphics]{par}} command to increase the margin size 
#'   (see the "mar" option in the \code{par} help page).
#'   
#'   The size of text labels can be modified by increasing or decreasing the
#'   \code{cex.main}, \code{cex.lab}, and \code{cex.axis} arguments:
#'   \itemize{
#'    \item{\code{cex.main}: }{controls the size of the plot title (specified 
#'                            in the \code{main} argument).}
#'    \item{\code{cex.lab}: }{controls the size of the axis labels on the
#'                           \emph{weighted degree}, \emph{node contribution},
#'                           and \emph{module summary} bar plots as well as
#'                           the size of the module labels and the heatmap 
#'                           legend titles.}
#'    \item{\code{cex.axis}: }{contols the size of the axis tick labels, 
#'                            including the node and sample labels.}
#'   }
#'   
#'   The position of these labels can be changed through the following 
#'   arguments:
#'   \itemize{
#'    \item{\code{xaxt.line}: }{controls the distance from the plot the x-axis
#'          tick labels are drawn on the \emph{module summary} bar plot.}
#'    \item{\code{xlab.line}: }{controls the distance from the plot the x-axis 
#'          label is drawn on the \emph{module summary} bar plot.}
#'    \item{\code{yaxt.line}: }{controls the distance from the plot the y-axis 
#'          tick labels are drawn on the \emph{weighted degree} and  
#'          \emph{node contribution} bar plots.}
#'    \item{\code{ylab.line}: }{controls the distance from the plot the y-axis
#'          label is drawn on the \emph{weighted degree} and 
#'          \emph{node contribution} bar plots.}
#'    \item{\code{main.line}: }{controls the distance from the plot the title is
#'          drawn.}
#'    \item{\code{naxt.line}: }{controls the distance from the plot the node 
#'          labels are drawn.}
#'    \item{\code{saxt.line}: }{controls the distance from the plot the sample 
#'          labels are drawn.}
#'    \item{\code{maxt.line}: }{controls the distance from the plot the module
#'          labels are drawn.}
#'    \item{\code{laxt.line}: }{controls the distance from the heatmap legends
#'          that the gradient legend labels are drawn.}
#'    \item{\code{legend.main.line}: }{controls the distance from the heatmap
#'          legends that the legend title is drawn.}
#'   }
#'   
#'   The rendering of node, sample, and module names can be disabled by setting
#'   \code{plotNodeNames}, \code{plotSampleNames}, and \code{plotModuleNames}
#'   to \code{FALSE}.
#'   
#'   The size of the axis ticks can be changed by increasing or decreasing the
#'   following arguments:
#'   \itemize{
#'    \item{\code{xaxt.tck}: }{size of the x-axis tick labels as a multiple of
#'          the height of the \emph{module summary} bar plot}
#'    \item{\code{yaxt.tck}: }{size of the y-axis tick labels as a multiple of 
#'          the width of the \emph{weighted degree} or \emph{node contribution}
#'          bar plots.}
#'    \item{\code{laxt.tck}: }{size of the heatmap legend axis ticks as a 
#'          multiple of the width of the data, correlation structure, or 
#'          network edge weight heatmaps.}
#'   }
#'   
#'   The placement of heatmap legends is controlled by the following arguments:
#'   \itemize{
#'    \item{\code{plotLegend}: }{if \code{FALSE} legend will not be drawn.}
#'    \item{\code{legend.position}: }{a multiple of the plot width, controls 
#'      the horizontal distance from the plot the legend is drawn.}
#'   } 
#'   
#'   The \code{drawBorders} argument controls whether borders are drawn around 
#'   the weighted degree, node contribution, or module summary bar plots. The 
#'   \code{lwd} argument controls the thickness of these borders, as well as 
#'   the thickness of axes and axis ticks.
#' }
#' \subsection{Modifying the color palettes:}{
#'   The \code{dataCols} and \code{dataRange} arguments control the appearance 
#'   of the data heatmap (see \code{\link{plotData}}). The gradient of colors 
#'   used on the heatmap can be changed by specifying a vector of colors to 
#'   interpolate between in \code{dataCols} and \code{dataRange} specifies the 
#'   range of values that maps to this gradient. Values outside of the 
#'   specified \code{dataRange} will be rendered with the colors used at either
#'   extreme of the gradient. The default gradient is determined based on the 
#'   \code{data} shown on the plot. If all values in the \code{data} matrix are
#'   positive, then the gradient is interpolated between white and green, where
#'   white is used for the smallest value and green for the largest. If all
#'   values are negative, then the gradient is interpolated between purple and
#'   white, where purple is used for the smallest value and white for the value
#'   closest to zero. If the data contains both positive and negative values, 
#'   then the gradient is interpolated between purple, white, and green, where 
#'   white is used for values of zero. In this case the range shown is always 
#'   centered at zero, with the values at either extreme determined by the 
#'   value in the rendered \code{data} with the strongest magnitude (the 
#'   maximum of the absolute value).
#'   
#'   The \code{corCols} and \code{corRange} arguments control the appearance of
#'   the correlation structure heatmap (see \code{\link{plotCorrelation}}). The
#'   gradient of colors used on the heatmap can be changed by specifying a
#'   vector of colors to interpolate between in \code{corCols}. By default,
#'   strong negative correlations are shown in blue, and strong positive
#'   correlations in red, and weak correlations as white. \code{corRange} 
#'   controls the range of values that this gradient maps to, by default, -1 to
#'   1. Changing this may be useful for showing differences where range of 
#'   correlation coefficients is small.
#'   
#'   The \code{netCols} and \code{netRange} arguments control the appearance of
#'   the network edge weight heatmap (see \code{\link{plotNetwork}}). The
#'   gradient of colors used on the heatmap can be changed by specifying a
#'   vector of colors to interpolate between in \code{netCols}. By default,
#'   weak or non-edges are shown in white, while strong edges are shown in red.
#'   The \code{netRange} controls the range of values this gradient maps to, 
#'   by default, 0 to 1. If \code{netRange} is set to \code{NA}, then the 
#'   gradient will be mapped to values between 0 and the maximum edge weight of
#'   the shown network.
#'   
#'   The \code{degreeCol} argument controls the color of the weighted degree
#'   bar plot (see \code{\link{plotDegree}}).
#'   
#'   The \code{contribCols} argument controls the color of the node 
#'   contribution bar plot (see \code{\link{plotContribution}}. This can be 
#'   specified as single value to be used for all nodes, or as two colors: one
#'   to use for nodes with positive contributions and one to use for nodes with
#'   negative contributions.
#'   
#'   The \code{summaryCols} argument controls the color of the module summary 
#'   bar plot (see \code{\link{plotSummary}}. This can be specified as single
#'   value to be used for all samples, or as two colors: one to use for samples
#'   with a positive module summary value and one fpr samples with a negative
#'   module summary value.
#'   
#'   The \code{naCol} argument controls the color of missing nodes and samples
#'   on the data, correlaton structure, and network edge weight heatmaps.
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
#' # Alternatively nodes can be ordered on the plot by degree in the test dataset
#' plotDegree(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4), discovery="discovery",
#'   test="test", orderNodesBy="test"
#' )
#' 
#' # Or by averaging the degree across datasets for a more robust ordering  
#' plotDegree(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4), discovery="discovery",
#'   test="test", orderNodesBy=c("discovery", "test")
#' )
#' 
#' # Arbitrary subsets can be plotted:
#' plotContribution(
#'   data=data_list[[1]][, 1:10], correlation=correlation_list[[1]][1:10, 1:10], 
#'   network=network_list[[1]][1:10, 1:10], orderNodesBy=NA
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
  backgroundLabel="0", discovery=NULL, test=NULL, verbose=TRUE,
  orderSamplesBy=NULL, orderNodesBy=NULL, orderModules=TRUE, plotNodeNames=TRUE, 
  plotSampleNames=TRUE, plotModuleNames=NULL, main="", main.line=1, lwd=1, 
  plotLegend=TRUE, legend.main="Data", legend.main.line=1, naxt.line=-0.5, 
  saxt.line=-0.5, maxt.line=3, legend.position=0.15, laxt.tck=0.03, laxt.line=3, 
  cex.axis=0.8, cex.lab=1.2, cex.main=2, dataCols=NULL, dataRange=NULL, 
  naCol="#bdbdbd", dryRun=FALSE
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters to catch errors prior to computation
  #-----------------------------------------------------------------------------
  
  old.par <- par(c("cex.axis", "cex.lab", "cex.main", "mar", "oma"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
    par(mar=old.par[[4]])
    par(oma=old.par[[5]])
    try(layout(1))
  }, add=TRUE)
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  if (is.null(data))
    stop("Cannot plot data matrix without 'data'")
  
  # Check plot-specific arguments
  checkPlotArgs(orderModules=orderModules, plotNodeNames=plotNodeNames, 
    plotSampleNames=plotSampleNames, plotModuleNames=plotModuleNames, 
    main=main, lwd=lwd, naxt.line=naxt.line, main.line=main.line, 
    saxt.line=saxt.line, maxt.line=maxt.line, laxt.line=laxt.line, 
    laxt.tck=laxt.tck,plotLegend=plotLegend, dataCols=dataCols,
    legend.main=legend.main, naCol=naCol, legend.position=legend.position, 
    dataRange=dataRange, dryRun=dryRun, legend.main.line=legend.main.line)
  
  # Handle variants that will not work for this plot function
  if (is.null(laxt.tck))
    stop("'laxt.tck' must be a numeric vector of length 1 or 'NA'")
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, plotFunction=TRUE, orderNodesBy, 
                         orderSamplesBy, orderModules)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  orderNodesBy <- finput$orderNodesBy
  orderSamplesBy <- finput$orderSamplesBy
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # Convert dataset indices to dataset names
  if (is.numeric(di))
    di <- datasetNames[di]
  if (is.numeric(ti))
    ti <- datasetNames[ti]

  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes, modules, and samples, identifying missing nodes and
  # samples on the plot, and get the network properties to be shown on the plot.
  #-----------------------------------------------------------------------------
  
  plotProps <- plotProps(data, correlation, network, moduleAssignments,
    modules, di, ti, orderNodesBy, orderSamplesBy, orderModules, datasetNames, 
    nDatasets, dryRun, verbose)
  testProps <- plotProps$testProps
  nodeOrder <- plotProps$nodeOrder
  moduleOrder <- plotProps$moduleOrder
  sampleOrder <- plotProps$sampleOrder
  na.pos.x <- plotProps$na.pos.x
  na.pos.y <- plotProps$na.pos.y
  presentNodes <- plotProps$presentNodes
  presentSamples <- plotProps$presentSamples
  
  #-----------------------------------------------------------------------------
  # Set default values for 'NULL' arguments
  #-----------------------------------------------------------------------------
  
  # Plot module names only if drawing more than one module
  if (is.null(plotModuleNames)) {
    plotModuleNames <- length(mods) > 1
  }
  
  # Set default color palettes for the data heatmap
  if (dryRun) {
    dat <- matrix(0, nrow=length(presentSamples), ncol=length(presentNodes))
    dimnames(dat) <- list(presentSamples, presentNodes)
    if (is.null(dataRange)) {
      dataRange <- c(-1, 1)
    }
  } else {
    dat <- data[[ti]][presentSamples, presentNodes] # also used for actual plot
    if (is.null(dataRange)) {
      dataRange <- range(dat)
      # Make sure the gradient is balanced around 0 if the default colors are
      # requested
      if (is.null(dataCols) && dataRange[1] < 0 && dataRange[2] > 0) {
        dataRange <- c(-1*max(abs(dataRange)), max(abs(dataRange)))
      }
    }
  }
  if (is.null(dataCols)) {
    if (all(dataRange >= 0)) {
      dataCols <- c("#FFFFFF", "#1B7837")
    } else if (all(dataRange <= 0)) {
      dataCols <- c("#762A83", "#FFFFFF")
    } else {
      dataCols <- c("#762A83", "#FFFFFF", "#1B7837")
    }
  }
  
  #-----------------------------------------------------------------------------
  # Plot the data matrix
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
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
    dat, dataCols, vlim=dataRange, legend.main.line=legend.main.line,
    moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.y, xaxt=xaxt, 
    yaxt=yaxt, plotLegend=plotLegend, main=main, main.line=main.line,
    legend.main=legend.main, plotModuleNames=plotModuleNames, 
    xaxt.line=naxt.line, yaxt.line=saxt.line, laxt.tck=laxt.tck,
    laxt.line=laxt.line, legend.line=legend.position, maxt.line=maxt.line,
    lwd=lwd, na.col=naCol, dryRun=dryRun
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
  backgroundLabel="0", discovery=NULL, test=NULL, verbose=TRUE,
  orderNodesBy=NULL, symmetric=FALSE, orderModules=TRUE, plotNodeNames=TRUE, 
  plotModuleNames=NULL, main="", main.line=1, lwd=1, plotLegend=TRUE, 
  legend.main="Correlation", legend.main.line=1, naxt.line=-0.5, maxt.line=3, 
  legend.position=NULL, laxt.tck=NULL, laxt.line=NULL, cex.axis=0.8, 
  cex.lab=1.2, cex.main=2, corCols=correlation.palette(), corRange=c(-1,1), 
  naCol="#bdbdbd", dryRun=FALSE
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters to catch errors prior to computation
  #-----------------------------------------------------------------------------
  
  old.par <- par(c("cex.axis", "cex.lab", "cex.main", "mar", "oma"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
    par(mar=old.par[[4]])
    par(oma=old.par[[5]])
    try(layout(1))
  }, add=TRUE)
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  # Check plot-specific arguments
  checkPlotArgs(orderModules=orderModules, plotNodeNames=plotNodeNames, 
    plotModuleNames=plotModuleNames, main=main, lwd=lwd, 
    naxt.line=naxt.line, maxt.line=maxt.line, laxt.line=laxt.line, 
    laxt.tck=laxt.tck, plotLegend=plotLegend, naCol=naCol, main.line=main.line,
    legend.main=legend.main, legend.position=legend.position, corCols=corCols, 
    corRange=corRange, symmetric=symmetric, dryRun=dryRun,
    legend.main.line=legend.main.line)

  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, plotFunction=TRUE, orderNodesBy, 
                         orderSamplesBy=NA, orderModules)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  orderNodesBy <- finput$orderNodesBy
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # Convert dataset indices to dataset names
  if (is.numeric(di))
    di <- datasetNames[di]
  if (is.numeric(ti))
    ti <- datasetNames[ti]
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and modules, identify missing nodes on the plot, and 
  # get the network properties to be shown on the plot.
  #-----------------------------------------------------------------------------
  
  plotProps <- plotProps(data, correlation, network, moduleAssignments,
    modules, di, ti, orderNodesBy, orderSamplesBy=NULL, orderModules, 
    datasetNames, nDatasets, dryRun, verbose)
  testProps <- plotProps$testProps
  nodeOrder <- plotProps$nodeOrder
  moduleOrder <- plotProps$moduleOrder
  na.pos.x <- plotProps$na.pos.x
  presentNodes <- plotProps$presentNodes
  
  #-----------------------------------------------------------------------------
  # Set default values for 'NULL' arguments
  #-----------------------------------------------------------------------------
  
  # Plot module names only if drawing more than one module
  if (is.null(plotModuleNames)) {
    plotModuleNames <- length(mods) > 1
  }
  
  # Default legend parameters depend on whether the heatmap is square or
  # triangular
  if (symmetric) {
    if (is.null(legend.position))
      legend.position <- 0.15
    if (is.null(laxt.line))
      laxt.line <- 3
    if (is.null(laxt.tck))
      laxt.tck <- 0.03
  } else {
    if (is.null(legend.position))
      legend.position <- 0.1
    if (is.null(laxt.line))
      laxt.line <- 2.5
    if (is.null(laxt.tck))
      laxt.tck <- 0.025
  }
  
  #-----------------------------------------------------------------------------
  # Plot the correlation structure
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  # Plot axis tick labels?
  naxt <- NULL
  if (plotNodeNames)
    naxt <- nodeOrder
  
  if (symmetric) {
    # Add space for the legend
    if (plotLegend) {
      old.mar <- par("mar")
      par(mar=par("mar")+c(0,0,0,3))
      on.exit({
        par(mar=old.mar)
      }, add=TRUE)
    }
    
    plotSquareHeatmap(
      correlation[[ti]][presentNodes, presentNodes], corCols, corRange, 
      moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.x, 
      xaxt=naxt, yaxt=naxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=naxt.line, yaxt.line=naxt.line, lwd=lwd,
      laxt.tck=laxt.tck, laxt.line=laxt.line, main.line=main.line,
      legend.line=legend.position, maxt.line=maxt.line, na.col=naCol,
      dryRun=dryRun,legend.main.line=legend.main.line
    )
  } else {
    plotTriangleHeatmap(
      correlation[[ti]][presentNodes, presentNodes], corCols, corRange, 
      moduleAssignments[[di]][nodeOrder], na.pos.x, xaxt=naxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=naxt.line,
      laxt.tck=laxt.tck, laxt.line=laxt.line, main.line=main.line,
      legend.line=legend.position, maxt.line=maxt.line, 
      lwd=lwd, na.col=naCol, dryRun=dryRun, legend.main.line=legend.main.line
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
  backgroundLabel="0", discovery=NULL, test=NULL, verbose=TRUE,
  orderNodesBy=NULL, symmetric=FALSE, orderModules=TRUE, plotNodeNames=TRUE, 
  plotModuleNames=NULL, main="", main.line=1, lwd=1, plotLegend=TRUE, 
  legend.main="Edge weight", legend.main.line=1, naxt.line=-0.5, maxt.line=3, 
  legend.position=NULL, laxt.tck=NULL, laxt.line=NULL, cex.axis=0.8, 
  cex.lab=1.2, cex.main=2, netCols=network.palette(), netRange=c(0,1), 
  naCol="#bdbdbd", dryRun=FALSE
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters to catch errors prior to computation
  #-----------------------------------------------------------------------------
  
  old.par <- par(c("cex.axis", "cex.lab", "cex.main", "mar", "oma"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
    par(mar=old.par[[4]])
    par(oma=old.par[[5]])
    try(layout(1))
  }, add=TRUE)
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  # Check plot-specific arguments
  checkPlotArgs(orderModules=orderModules, plotNodeNames=plotNodeNames, 
    plotModuleNames=plotModuleNames, main=main, lwd=lwd, 
    naxt.line=naxt.line, maxt.line=maxt.line, laxt.line=laxt.line, 
    laxt.tck=laxt.tck, plotLegend=plotLegend, naCol=naCol, main.line=main.line,
    legend.main=legend.main, legend.position=legend.position,
    symmetric=symmetric, netCols=netCols, netRange=netRange, dryRun=dryRun,
    legend.main.line=legend.main.line)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, plotFunction=TRUE, orderNodesBy, 
                         orderSamplesBy=NA, orderModules)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  orderNodesBy <- finput$orderNodesBy
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # Convert dataset indices to dataset names
  if (is.numeric(di))
    di <- datasetNames[di]
  if (is.numeric(ti))
    ti <- datasetNames[ti]
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and modules, identify missing nodes on the plot, and 
  # get the network properties to be shown on the plot.
  #-----------------------------------------------------------------------------
  
  plotProps <- plotProps(data, correlation, network, moduleAssignments,
    modules, di, ti, orderNodesBy, orderSamplesBy=NULL, orderModules, 
    datasetNames, nDatasets, dryRun, verbose)
  testProps <- plotProps$testProps
  nodeOrder <- plotProps$nodeOrder
  moduleOrder <- plotProps$moduleOrder
  na.pos.x <- plotProps$na.pos.x
  presentNodes <- plotProps$presentNodes
  
  #-----------------------------------------------------------------------------
  # Set default values for 'NULL' arguments
  #-----------------------------------------------------------------------------
  
  # Plot module names only if drawing more than one module
  if (is.null(plotModuleNames)) {
    plotModuleNames <- length(mods) > 1
  }
  
  # Default legend parameters depend on whether the heatmap is square or
  # triangular
  if (symmetric) {
    if (is.null(legend.position))
      legend.position <- 0.15
    if (is.null(laxt.line))
      laxt.line <- 3
    if (is.null(laxt.tck))
      laxt.tck <- 0.03
  } else {
    if (is.null(legend.position))
      legend.position <- 0.1
    if (is.null(laxt.line))
      laxt.line <- 2.5
    if (is.null(laxt.tck))
      laxt.tck <- 0.025
  }
  
  #-----------------------------------------------------------------------------
  # Plot the network edge weights
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  # Plot axis tick labels?
  naxt <- NULL
  if (plotNodeNames)
    naxt <- nodeOrder
  
  if (symmetric) {
    # Add space for the legend
    if (plotLegend) {
      old.mar <- par("mar")
      par(mar=par("mar")+c(0,0,0,3))
      on.exit({
        par(mar=old.mar)
      }, add=TRUE)
    }
    
    plotSquareHeatmap(
      network[[ti]][presentNodes, presentNodes], netCols, netRange, 
      moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.x, 
      xaxt=naxt, yaxt=naxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=naxt.line, yaxt.line=naxt.line, 
      laxt.tck=laxt.tck, laxt.line=laxt.line, main.line=main.line,
      legend.line=legend.position, maxt.line=maxt.line,
      lwd=lwd, na.col=naCol, dryRun=dryRun, legend.main.line=legend.main.line
    )
  } else {
    plotTriangleHeatmap(
      network[[ti]][presentNodes, presentNodes], netCols, netRange, 
      moduleAssignments[[di]][nodeOrder], na.pos.x,xaxt=naxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=naxt.line,
      laxt.tck=laxt.tck, laxt.line=laxt.line, main.line=main.line, 
      legend.line=legend.position, maxt.line=maxt.line, 
      lwd=lwd, na.col=naCol, dryRun=dryRun, legend.main.line=legend.main.line
    )
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotContribution}: Plot a bar chart of the node contribution for
#' one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotContribution <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, verbose=TRUE,
  orderNodesBy=NULL, orderModules=TRUE, plotNodeNames=TRUE, 
  plotModuleNames=NULL, main="", main.line=1, ylab.line=2.5, lwd=1, 
  drawBorders=FALSE, naxt.line=-0.5, maxt.line=3, yaxt.line=0, yaxt.tck=-0.035, 
  cex.axis=0.8, cex.lab=1.2, cex.main=2, contribCols=c("#A50026", "#313695"), 
  naCol="#bdbdbd", dryRun=FALSE  
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters to catch errors prior to computation
  #-----------------------------------------------------------------------------
  
  old.par <- par(c("cex.axis", "cex.lab", "cex.main", "mar", "oma"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
    par(mar=old.par[[4]])
    par(oma=old.par[[5]])
    try(layout(1))
  }, add=TRUE)
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  if (is.null(data))
    stop("Cannot plot node contribution without 'data'")
  
  # Check plot-specific arguments
  checkPlotArgs(orderModules=orderModules, plotNodeNames=plotNodeNames, 
    plotModuleNames=plotModuleNames, main=main, lwd=lwd, main.line=main.line,
    ylab.line=ylab.line, drawBorders=drawBorders, naxt.line=naxt.line, 
    maxt.line=maxt.line, contribCols=contribCols, naCol=naCol, dryRun=dryRun, 
    yaxt.line=yaxt.line, yaxt.tck=yaxt.tck)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, plotFunction=TRUE, orderNodesBy, 
                         orderSamplesBy=NA, orderModules)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  orderNodesBy <- finput$orderNodesBy
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # Convert dataset indices to dataset names
  if (is.numeric(di))
    di <- datasetNames[di]
  if (is.numeric(ti))
    ti <- datasetNames[ti]
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and modules, identify missing nodes on the plot, and 
  # get the network properties to be shown on the plot.
  #-----------------------------------------------------------------------------
  
  plotProps <- plotProps(data, correlation, network, moduleAssignments,
    modules, di, ti, orderNodesBy, orderSamplesBy=NULL, orderModules, 
    datasetNames, nDatasets, dryRun, verbose)
  testProps <- plotProps$testProps
  nodeOrder <- plotProps$nodeOrder
  moduleOrder <- plotProps$moduleOrder
  na.pos.x <- plotProps$na.pos.x
  presentNodes <- plotProps$presentNodes
  
  #-----------------------------------------------------------------------------
  # Set default values for 'NULL' arguments
  #-----------------------------------------------------------------------------
  
  # Plot module names only if drawing more than one module
  if (is.null(plotModuleNames)) {
    plotModuleNames <- length(mods) > 1
  }

  #-----------------------------------------------------------------------------
  # Plot the Node contribution
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  if (dryRun) {
    nodeContribVec <- rep(0, length(nodeOrder))
    names(nodeContribVec) <- nodeOrder
  } else {
    nodeContribVec <- foreach(mi = moduleOrder, .combine=c) %do% {
      testProps[[mi]]$contribution
    }
    nodeContribVec <- nodeContribVec[nodeOrder]
  }
  
  # Plot bar chart
  plotBar(
    nodeContribVec, c(-1,1), moduleAssignments[[di]][nodeOrder],
    contribCols, drawBorders=drawBorders, main.line=main.line, 
    xaxt=plotNodeNames, plotModuleNames=plotModuleNames, ylab.line=ylab.line,
    xaxt.line=naxt.line, maxt.line=maxt.line, yaxt.line=yaxt.line, 
    yaxt.tck=yaxt.tck, main=main, ylab="Node contribution", 
    lwd=lwd, na.col=naCol, dryRun=dryRun
  )
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}

#' \code{plotDegree:} Plot a bar chart of the scaled weighted degree 
#' (see details) for one or more network modules.
#' 
#' @rdname plotTopology
#' @export
plotDegree <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, verbose=TRUE,
  orderNodesBy=NULL, orderModules=TRUE, plotNodeNames=TRUE, 
  plotModuleNames=NULL, main="", main.line=1, lwd=1, drawBorders=FALSE, 
  naxt.line=-0.5, maxt.line=3, yaxt.line=0, yaxt.tck=-0.035, ylab.line=2.5, 
  cex.axis=0.8, cex.lab=1.2, cex.main=2, degreeCol="#feb24c", naCol="#bdbdbd", 
  dryRun=FALSE
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters to catch errors prior to computation
  #-----------------------------------------------------------------------------
  
  old.par <- par(c("cex.axis", "cex.lab", "cex.main", "mar", "oma"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
    par(mar=old.par[[4]])
    par(oma=old.par[[5]])
    try(layout(1))
  }, add=TRUE)
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  # Check plot-specific arguments
  checkPlotArgs(orderModules=orderModules, plotNodeNames=plotNodeNames, 
    plotModuleNames=plotModuleNames, main=main, lwd=lwd, main.line=main.line,
    drawBorders=drawBorders, naxt.line=naxt.line, maxt.line=maxt.line, 
    degreeCol=degreeCol, naCol=naCol, dryRun=dryRun, yaxt.line=yaxt.line, 
    yaxt.tck=yaxt.tck, ylab.line=ylab.line)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, plotFunction=TRUE, orderNodesBy, 
                         orderSamplesBy=NA, orderModules)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  orderNodesBy <- finput$orderNodesBy
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # Convert dataset indices to dataset names
  if (is.numeric(di))
    di <- datasetNames[di]
  if (is.numeric(ti))
    ti <- datasetNames[ti]

  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes and modules, identify missing nodes on the plot, and 
  # get the network properties to be shown on the plot.
  #-----------------------------------------------------------------------------
  
  plotProps <- plotProps(data, correlation, network, moduleAssignments,
    modules, di, ti, orderNodesBy, orderSamplesBy=NULL, orderModules, 
    datasetNames, nDatasets, dryRun, verbose)
  testProps <- plotProps$testProps
  nodeOrder <- plotProps$nodeOrder
  moduleOrder <- plotProps$moduleOrder
  na.pos.x <- plotProps$na.pos.x
  presentNodes <- plotProps$presentNodes
  
  #-----------------------------------------------------------------------------
  # Set default values for 'NULL' arguments
  #-----------------------------------------------------------------------------
  
  # Plot module names only if drawing more than one module
  if (is.null(plotModuleNames)) {
    plotModuleNames <- length(mods) > 1
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Node degree
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  if (dryRun) {
    wDegreeVec <- rep(0, length(nodeOrder))
    names(wDegreeVec) <- nodeOrder
  } else {
    # (Normalised) weighted degree vector
    wDegreeVec <- foreach(mi = moduleOrder, .combine=c) %do% {
      testProps[[mi]]$degree/max(na.omit(testProps[[mi]]$degree))
    }
    wDegreeVec <- wDegreeVec[nodeOrder]
  }
  
  # Plot bar chart
  plotBar(
    wDegreeVec, c(0,1), moduleAssignments[[di]][nodeOrder],
    degreeCol, drawBorders=drawBorders, main.line=main.line, 
    xaxt=plotNodeNames, plotModuleNames=plotModuleNames, 
    xaxt.line=naxt.line, maxt.line=maxt.line, ylab.line=ylab.line,
    yaxt.line=yaxt.line, yaxt.tck=yaxt.tck, main=main,
    ylab="Weighted degree", lwd=lwd, na.col=naCol,
    dryRun=dryRun
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
  backgroundLabel="0", discovery=NULL, test=NULL, verbose=TRUE,
  orderSamplesBy=NULL, orderNodesBy=NULL, orderModules=TRUE, 
  plotSampleNames=TRUE, plotModuleNames=NULL, main="", main.line=1, 
  xlab.line=2.5, lwd=1, drawBorders=FALSE, saxt.line=-0.5, maxt.line=0, 
  xaxt.line=0, xaxt.tck=-0.025, cex.axis=0.8, cex.lab=1.2, cex.main=2, 
  summaryCols=c("#1B7837", "#762A83"), naCol="#bdbdbd", dryRun=FALSE
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters to catch errors prior to computation
  #-----------------------------------------------------------------------------
  
  old.par <- par(c("cex.axis", "cex.lab", "cex.main", "mar", "oma"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
    par(mar=old.par[[4]])
    par(oma=old.par[[5]])
    try(layout(1))
  }, add=TRUE)
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  if (is.null(data))
    stop("Cannot plot module summary without 'data'")
  
  # Check plot-specific arguments
  checkPlotArgs(orderModules=orderModules, plotSampleNames=plotSampleNames, 
    plotModuleNames=plotModuleNames, main=main, lwd=lwd, main.line=main.line,
    saxt.line=saxt.line, maxt.line=maxt.line, summaryCols=summaryCols, 
    naCol=naCol, dryRun=dryRun, xaxt.line=xaxt.line, xaxt.tck=xaxt.tck,
    xlab.line=xlab.line)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, plotFunction=TRUE, orderNodesBy, 
                         orderSamplesBy, orderModules)
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  correlation <- finput$correlation
  network <- finput$network
  moduleAssignments <- finput$moduleAssignments
  modules <- finput$modules
  nDatasets <- finput$nDatasets
  datasetNames <- finput$datasetNames
  orderSamplesBy <- finput$orderSamplesBy
  orderNodesBy <- finput$orderNodesBy
  
  # Indexes for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE
  
  # Convert dataset indices to dataset names
  if (is.numeric(di))
    di <- datasetNames[di]
  if (is.numeric(ti))
    ti <- datasetNames[ti]
  
  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes, modules, and samples, identifying missing nodes and
  # samples on the plot, and get the network properties to be shown on the plot.
  #-----------------------------------------------------------------------------
  
  plotProps <- plotProps(data, correlation, network, moduleAssignments,
    modules, di, ti, orderNodesBy, orderSamplesBy, orderModules, datasetNames, 
    nDatasets, dryRun, verbose)
  testProps <- plotProps$testProps
  nodeOrder <- plotProps$nodeOrder
  moduleOrder <- plotProps$moduleOrder
  sampleOrder <- plotProps$sampleOrder
  na.pos.x <- plotProps$na.pos.x
  na.pos.y <- plotProps$na.pos.y
  presentNodes <- plotProps$presentNodes
  presentSamples <- plotProps$presentSamples
  
  #-----------------------------------------------------------------------------
  # Set default values for 'NULL' arguments
  #-----------------------------------------------------------------------------
  
  # Plot module names only if drawing more than one module
  if (is.null(plotModuleNames)) {
    plotModuleNames <- length(mods) > 1
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Node contribution
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  
  if (dryRun) {
    summaries=matrix(0, ncol=length(mods), nrow=length(sampleOrder))
    dimnames(summaries) <- list(sampleOrder, mods)
    summaries.range <- rep(list(c(-1, 1)), length(mods))
    names(summaries.range) <- mods
  } else {
    # Summary profile matrix
    summaries <- foreach(mi = moduleOrder, .combine=cbind) %do% {
      matrix(
        insert.nas(testProps[[mi]]$summary[presentSamples], na.pos.y),
        ncol=1
      )
    }
    colnames(summaries) <- moduleOrder
    rownames(summaries) <- sampleOrder
    summaries.range <- lapply(1:ncol(summaries), function(ii) { 
      range(summaries[,ii])
    })
  }
  
  # Plot bar chart
  plotMultiBar(
    summaries, summaries.range, main.line=main.line, xlab.line=xlab.line,
    cols=summaryCols, drawBorders=drawBorders, main=main, 
    yaxt=plotSampleNames, plotModuleNames=plotModuleNames, yaxt.line=saxt.line, 
    xaxt.line=xaxt.line, xaxt.tck=xaxt.tck, maxt.line=maxt.line, 
    xlab="Module summary", lwd=lwd, na.col=naCol, 
    dryRun=dryRun
  )
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}
