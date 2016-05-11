#' Plot the topology of a network module
#' 
#' @description 
#' Plot the correlation structure, network edges, scaled weighted degree, 
#' node contribtuion, module data, and module summary vectors of one or
#' more network modules.
#' 
#' Individual components of the module plot can be plotted using 
#' \code{\link{plotCorrelation}}, \code{\link{plotNetwork}}, 
#' \code{\link{plotDegree}}, \code{\link{plotContribution}}, 
#' \code{\link{plotData}}, and \code{\link{plotSummary}}.
#' 
#' @inheritParams common_params
#' @inheritParams par_param
#' @inheritParams orderModules_param
#' @inheritParams plot_params
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
#'     \item{\code{network}: }{
#'       a list of interaction networks, one for each dataset.
#'     }
#'     \item{\code{data}: }{
#'       a list of data matrices used to infer those networks, one for each 
#'       dataset.
#'     }
#'     \item{\code{correlation}: }{
#'      a list of matrices containing the pairwise correlation coefficients 
#'      between variables/nodes in each dataset.
#'     } 
#'     \item{\code{moduleAssignments}: }{
#'      a list of vectors, one for each \emph{discovery} dataset, containing 
#'      the module assignments for each node in that dataset.
#'     }
#'     \item{\code{modules}: }{
#'      a list of vectors, one for each \emph{discovery} dataset, containing
#'      the names of the modules from that dataset to analyse.  
#'     }
#'     \item{\code{discovery}: }{
#'       a vector indicating the names or indices of the previous arguments' 
#'       lists to use as the \emph{discovery} dataset(s) for the analyses.
#'     }
#'     \item{\code{test}: }{
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
#'   By default, samples in the data heatmap and accompanying module summary bar
#'   plot are ordered in descending order of \emph{module summary} in the drawn 
#'   dataset (specified by the \code{test} argument). If multiple modules are 
#'   drawn, samples are ordered as per the left-most module on the plot.
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
#'   For optimal results we recommend viewing single modules on a PNG device
#'   with a width of 1500, a height of 2700 and a nominal resolution of 300 
#'   (\code{png(filename, width=5*300, height=9*300, res=300))}).
#'   
#'   \strong{Warning}: PDF and other vectorized devices should not be used when
#'   plotting more than a hundred nodes. Large files will be generated which
#'   may cause image editing programs such as Inkscape or Illustrator to crash
#'   when polishing figures for publication.
#'   
#'   When \code{dryRun} is \code{TRUE} only the axes, legends, labels, and
#'   title will be drawn, allowing for quick iteration of customisable
#'   parameters to get the plot layout correct.
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
#' \subsection{Embedding in Rmarkdown documents}{
#'   The chunk option \code{fig.keep="last"} should be set to avoid an empty 
#'   plot being embedded above the plot generated by \code{plotModule}. This
#'   empty plot is generated so that an error will be thrown as early as 
#'   possible if the margins are too small to be displayed. Normally, these
#'   are drawn over with the actual plot components when drawing the plot on
#'   other graphical devices.
#' }
#' 
#' @seealso
#' \code{\link{plotCorrelation}}, 
#' \code{\link{plotNetwork}},
#' \code{\link{plotDegree}},
#' \code{\link{plotContribution}},
#' \code{\link{plotData}}, and
#' \code{\link{plotSummary}}.
#' 
#' @examples
#' \dontrun{
#' # load in example data, correlation, and network matrices for a discovery 
#' # and test dataset:
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
#' # Plot module 1, 2 and 4 in the discovery dataset
#' plotModule(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4)
#' )
#' 
#' # Now plot them in the test dataset (module 2 does not replicate)
#' plotModule(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 2, 4), discovery="discovery",
#'   test="test"
#' )
#' 
#' # Plot modules 1 and 4, which replicate, in the test datset ordering nodes
#' # by weighted degree averaged across the two datasets
#' plotModule(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules=c(1, 4), discovery="discovery",
#'   test="test", orderNodesBy=c("discovery", "test")
#' )
#' }
#' 
#' @name plotModule
#' @export
plotModule <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, verbose=TRUE,
  orderSamplesBy=NULL, orderNodesBy=NULL, orderModules=TRUE, plotNodeNames=TRUE, 
  plotSampleNames=TRUE, plotModuleNames=NULL, main="Module Topology", 
  main.line=1, drawBorders=FALSE, lwd=1, naxt.line=-0.5, saxt.line=-0.5, 
  maxt.line=NULL, xaxt.line=-0.5, xaxt.tck=-0.025, xlab.line=2.5, yaxt.line=0, 
  yaxt.tck=-0.15, ylab.line=2.5, laxt.line=2.5, laxt.tck=0.04, cex.axis=0.8, 
  legend.main.line=1.5, cex.lab=1.2, cex.main=2, dataCols=NULL, dataRange=NULL, 
  corCols=correlation.palette(), corRange=c(-1,1), netCols=network.palette(), 
  netRange=c(0,1), degreeCol="#feb24c", contribCols=c("#A50026", "#313695"), 
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
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  checkPlotArgs(orderModules=orderModules, plotNodeNames=plotNodeNames, 
    plotSampleNames=plotSampleNames, plotModuleNames=plotModuleNames, 
    main=main, drawBorders=drawBorders, lwd=lwd, naxt.line=naxt.line, 
    saxt.line=saxt.line, maxt.line=maxt.line, xaxt.line=xaxt.line, 
    xaxt.tck=xaxt.tck, yaxt.line=yaxt.line, yaxt.tck=yaxt.tck, 
    laxt.tck=laxt.tck, laxt.line=laxt.line, ylab.line=ylab.line, 
    xlab.line=xlab.line, main.line=main.line, dataCols=dataCols, 
    dataRange=dataRange, corCols=corCols, corRange=corRange, netCols=netCols, 
    netRange=netRange, degreeCol=degreeCol, contribCols=contribCols, 
    summaryCols=summaryCols, naCol=naCol, legend.main.line=legend.main.line, 
    dryRun=dryRun)
  
  # Handle variants that will not work for this plot function
  if (is.null(laxt.tck))
    stop("'laxt.tck' must be a numeric vector of length 1 or 'NA'")
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir, plotFunction=TRUE, orderNodesBy, 
                         orderSamplesBy, orderModules, dryRun)
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
  orderNodesBy <- finput$orderNodesBy
  orderSamplesBy <- finput$orderSamplesBy
  
  # Indices for this function
  di <- finput$discovery
  ti <- finput$test[[di]]
  mods <- modules[[di]]
  mi <- NULL # initialise to suppress CRAN NOTE

  # Convert dataset indices to dataset names
  if (is.numeric(di))
    di <- datasetNames[di]
  if (is.numeric(ti))
    ti <- datasetNames[ti]
  
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  # Create empty plot windows so that we fail quickly if the margins are too
  # large
  if (is.null(scaledData[[ti]])) {
    layout(mat=matrix(1:3, ncol=1), heights=c(0.4, 0.4, 0.2))
  } else {
    summary.window <- min(0.2 + (length(mods) - 1)* 0.1, 0.5) 
    layout(
      mat=matrix(c(rep(8, 5), 7, 1:6), ncol=2), 
      heights=c(0.7/3, 0.7/3, 0.12, 0.12, 0.06, 0.7/3),
      widths=c(summary.window, 1-summary.window)
    )
  }
  par(oma=old.par[["mar"]]+old.par[["oma"]])
  
  tryCatch({
    par(mar=c(1, 1, 1, 1))
    emptyPlot(xlim=c(0, 1), ylim=c(0, 1), bty="n") # Correlation heatmap plot
    par(mar=c(1, 1, 1, 1))
    emptyPlot(xlim=c(0, 1), ylim=c(0, 1), bty="n") # Network heatmap
    par(mar=c(1, 1, 1, 1))
    emptyPlot(xlim=c(0, 1), ylim=c(0, 1), bty="n") # Degree barplot
    if (!is.null(scaledData[[ti]])) {
      par(mar=c(1, 1, 1, 1))
      emptyPlot(xlim=c(0, 1), ylim=c(0, 1), bty="n") # Contribution barplot
      par(mar=c(1,1,1,1))
      emptyPlot(c(0,1), c(0,1), bty="n") # Empty space for data heatmap legend
      par(mar=c(1, 1, 1, 1))
      emptyPlot(c(0,1), c(0,1), bty="n") # Data heatmap
      par(mar=c(1, 1, 1, 1))
      emptyPlot(xlim=c(0, 1), ylim=c(0, 1), bty="n") # Summary bar plots
      par(mar=c(0, 0, 0, 0))
      emptyPlot(xlim=c(0, 1), ylim=c(0, 1), bty="n") # Empty region layout
    }
  }, error=function(e) {
    stop("graphics device too small to render plot")
  })

  vCat(verbose, 0, "User input ok!")
  
  #-----------------------------------------------------------------------------
  # Get ordering of nodes, modules, and samples, identifying missing nodes and
  # samples on the plot, and get the network properties to be shown on the plot.
  #-----------------------------------------------------------------------------
  
  plotProps <- plotProps(scaledData, correlation, network, moduleAssignments,
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
  # Get the data to be shown on the plot
  #-----------------------------------------------------------------------------
  
  if (dryRun) {
    # In this case these properties have not been calculated, so just make some
    # placehold vectors so that the plot functions work
    wDegreeVec <- rep(0, length(nodeOrder))
    names(wDegreeVec) <- nodeOrder
    if (!is.null(scaledData[[ti]])) {
      nodeContribVec <- rep(0, length(nodeOrder))
      names(nodeContribVec) <- nodeOrder
      
      summaries=matrix(0, ncol=length(mods), nrow=length(sampleOrder))
      dimnames(summaries) <- list(sampleOrder, mods)
      summaries.range <- rep(list(c(-1, 1)), length(mods))
      names(summaries.range) <- mods
    }
  } else {
    # (Normalised) weighted degree vector
    wDegreeVec <- foreach(mi = moduleOrder, .combine=c) %do% {
      testProps[[mi]]$degree/max(na.omit(testProps[[mi]]$degree))
    }
    wDegreeVec <- wDegreeVec[nodeOrder]
    
    if (!is.null(scaledData[[ti]])) {
      # node contribution
      nodeContribVec <- foreach(mi = moduleOrder, .combine=c) %do% {
        testProps[[mi]]$contribution
      }
      nodeContribVec <- nodeContribVec[nodeOrder]
      
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
  }
    
  #-----------------------------------------------------------------------------
  # Set default values for 'NULL' arguments
  #-----------------------------------------------------------------------------
  
  # Plot module names only if drawing more than one module
  if (is.null(plotModuleNames)) {
    plotModuleNames <- length(mods) > 1
  }
  
  # Set the location of module names in the margin based on whether or not 
  # we're plotting the node names
  if (is.null(maxt.line) && !plotNodeNames) {
    maxt.line <- -0.5
  } else if (is.null(maxt.line) && plotNodeNames) {
    maxt.line <- 3
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
  
  naxt <- NULL
  if (plotNodeNames)
    naxt <- nodeOrder
  
  #-----------------------------------------------------------------------------
  # Plot topology components
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "rendering plot components...")
  # Plot correlation
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    correlation[[ti]][presentNodes, presentNodes], 
    corCols, corRange, moduleAssignments[[di]][nodeOrder], na.pos.x, 
    plotLegend=TRUE, main=main, main.line=main.line, legend.main="Correlation", 
    plotModuleNames=FALSE, laxt.tck=laxt.tck, laxt.line=laxt.line, na.col=naCol,
    legend.line=0.1, maxt.line=maxt.line, lwd=lwd, 
    legend.main.line=legend.main.line, dryRun=dryRun
  )

  # Plot network
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    network[[ti]][presentNodes, presentNodes], netCols, netRange,
    moduleAssignments[[di]][nodeOrder], na.pos.x, plotLegend=TRUE, main="", 
    legend.main="Edge weights", plotModuleNames=FALSE, 
    laxt.tck=laxt.tck, na.col=naCol, legend.main.line=legend.main.line,
    laxt.line=laxt.line, legend.line=0.1, maxt.line=maxt.line,
    lwd=lwd, dryRun=dryRun
  )
  
  # Plot weighted degree
  par(mar=c(1,1,1,1))
  if (is.null(scaledData[[ti]])) {
    plotBar(
      wDegreeVec, c(0,1), moduleAssignments[[di]][nodeOrder], degreeCol, 
      drawBorders=drawBorders, plotModuleNames=plotModuleNames, 
      xaxt=plotNodeNames, xaxt.line=naxt.line, main="",
      ylab="Weighted\ndegree", maxt.line=maxt.line, yaxt.line=yaxt.line,
      yaxt.tck=yaxt.tck, lwd=lwd, na.col=naCol, dryRun=dryRun
    )
  } else {
    plotBar(
      wDegreeVec, c(0,1), moduleAssignments[[di]][nodeOrder], degreeCol, 
      drawBorders=drawBorders, plotModuleNames=FALSE, main="", xaxt=FALSE,
      ylab="Weighted\ndegree", maxt.line=maxt.line, yaxt.line=yaxt.line,
      yaxt.tck=yaxt.tck, ylab.line=ylab.line, na.col=naCol, lwd=lwd, 
      dryRun=dryRun
    )
  }
  
  if (!is.null(scaledData[[ti]])) {
    # Plot Module Membership
    par(mar=c(1, 1, 1, 1))
    plotBar(
      nodeContribVec, c(-1,1), moduleAssignments[[di]][nodeOrder], 
      contribCols, drawBorders=drawBorders, plotModuleNames=FALSE, main="", 
      xaxt=FALSE, ylab="Node\ncontribution", maxt.line=maxt.line, 
      yaxt.line=yaxt.line, yaxt.tck=yaxt.tck, ylab.line=ylab.line, na.col=naCol,
      lwd=lwd, dryRun=dryRun
    )
    
    # Plot the data matrix
    par(mar=c(1,1,1,1))
    emptyPlot(c(0,1), c(0,1), bty="n")

    yaxt <- NULL
    if (plotSampleNames)
      yaxt <- sampleOrder
    par(mar=c(1, 1, 1, 1))
    plotSquareHeatmap(
      dat, dataCols, vlim=dataRange,
      moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.y, 
      xaxt=naxt, yaxt=NULL, plotLegend=FALSE, main="",
      legend.main="", plotModuleNames=plotModuleNames,
      xaxt.line=naxt.line, maxt.line=maxt.line, lwd=lwd,
      dryRun=dryRun
    )
    
    # Plot data legend
    xHalfUnit <- (1/(ncol(dat) + length(na.pos.x)))/2
    yHalfUnit <- (1/(nrow(dat) + length(na.pos.y)))/2
    addGradientLegend(
      dataCols, dataRange, TRUE, main="Module data",
      xlim=c(xHalfUnit + 0.1, 1 - xHalfUnit - 0.1), 
      ylim=c(1 + yHalfUnit + 0.2, 1 + yHalfUnit + 0.3),  
      tck=laxt.tck, lwd=lwd, legend.main.line=legend.main.line,
      axis.line=laxt.line, srt=0
    )
    
    # Plot bar chart
    xlab <- "Module summary"
    if (length(mods) == 1) 
      xlab <- gsub(" ", "\n", xlab)
    par(mar=c(1, 1, 1, 1))
    plotMultiBar(
      summaries, summaries.range,
      cols=summaryCols , drawBorders=drawBorders, lwd=lwd,
      yaxt=plotSampleNames, plotModuleNames=plotModuleNames, 
      yaxt.line=saxt.line, maxt.line=0, xlab=xlab, xaxt.line=xaxt.line, 
      xlab.line=xlab.line, xaxt.tck=xaxt.tck, cex.modules=par("cex.lab"), 
      na.col=naCol, dryRun=dryRun
    )
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
}