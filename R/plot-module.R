#' Plot the topology of a network module
#' 
#' @description 
#' Plot the correlation structure, network edges, (normalised) connectivity, 
#' module membership, underlying data, and module summary vector of one or
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
#' @param main title for the plot.
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
#' \subsection{Plot customisation:}{
#'   Although reasonable default values for most parameters have been provided,
#'   the rendering of axes and titles may need adjusting depending on the size
#'   of the plot window. The parameters \code{gaxt.line}, \code{saxt.line}, 
#'   \code{maxt.line}, and \code{laxt.line} control the distance from each plot
#'   window that the node labels, sample labels, module labels, and legend 
#'   labels are rendered. 
#'   
#'   \code{legend.tick.size} controls the length of the 
#'   axis ticks on each of the legends relative to the correlation, edge weight,
#'   and data matrix heatmap plot windows. 
#'   
#'   \code{cex.main} controls the relative text size of the plot title
#'   (specified by the \code{main} argument). \code{cex.axis} controls the
#'   relative text size of the node and sample labels. \code{cex.lab} controls
#'   the relative text size of the bar plot axis labels, module labels, and the
#'   legend titles.
#'   
#'   The rendering of node, sample, and module names can be disabled by setting
#'   \code{plotNodeNames}, \code{plotSampleNames}, and \code{plotModuleNames} to
#'   \code{FALSE}.
#'   
#'   The \code{drawBorders} argument controls whether borders are drawn around
#'   the connectivity, module membership, or module summary bar plots.
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
#' plotModule(geA, coexpA, adjA, moduleAssignments, modules="2")
#' 
#' ## Example 2: Plot an arbitrary set of genes in cohort A
#' plotModule(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
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
#' data <- list(
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
#' plotModule(
#'   data, correlation, network, moduleAssignments,
#'   modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' }
#' 
#' @name plotModule
#' @export
plotModule <- function(
  data, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=1, test=1, nCores=1, verbose=TRUE,
  orderSamplesBy="test", orderNodesBy="discovery",
  orderModules=TRUE, plotNodeNames=TRUE, plotSampleNames=TRUE, plotModuleNames,
  main="Module Topology", drawBorders=FALSE, border.width=2, gaxt.line=-0.5, 
  saxt.line=-0.5, maxt.line, legend.tick.size=0.04, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  
  # Definition in here to prevent code potentially breaking elsewhere.
  is.vector <- function(obj) {
    base::is.vector(obj) & !is.list(obj)
  }
  
  #-----------------------------------------------------------------------------
  # Set graphical parameters
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
    layout(1)
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
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
  if (!is.logical(orderModules) | is.na(orderModules) | length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # At this time, we can only plot within one dataset.
  if (!is.vector(discovery) | !is.vector(test) | 
      length(discovery) > 1 | length(test) > 1) {
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
  
  # set up 'discovery' as 'test' so we can use it on 'netPropsInternal'
  discAsTest <- list(discovery)
  names(discAsTest) <- discovery

  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  if (missing(plotModuleNames))
    plotModuleNames <- length(mods) > 1
  
  if ((orderSamplesBy == "discovery" & is.null(scaledData[[di]])) |
      (orderSamplesBy == "test" & is.null(scaledData[[ti]]))) {
    stop("'data' not provided for 'orderSamplesBy' dataset") 
  }
  
  if (orderSamplesBy == "discovery" & 
      sum(rownames(scaledData[di]) %in% rownames(scaledData[[ti]])) == 0) {
    stop("'orderBySamples' can only be ", '"discovery"', " when the same",
         " samples are present in both the 'discovery' and 'test' datasets")
  }
  
  if ((orderModules & length(mods) > 1) & 
      ((orderNodesBy == "discovery" & is.null(scaledData[[di]])) | 
       (orderNodesBy == "test" & is.null(scaledData[[ti]])))) {
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
  if (orderNodesBy == "discovery" & di != ti) {
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
  if (orderSamplesBy == "discovery" & di != ti) {
    # This skips all of the data verification
    if (!exists("discProps")) {
      discProps <- netPropsInternal(
        scaledData, correlation, network, moduleAssignments, 
        modules, discovery, discAsTest, datasetNames, FALSE
      )
    }
    sampleOrder <- sampleOrderInternal(discProps, verbose, FALSE)[[di]][[di]][[1]]
  } 
  # Case 2: order samples as they're provided by the user
  else if (orderSamplesBy == "none") {
    sampleOrder <- seq_along(testProps[[di]][[ti]][[moduleOrder[1]]]$summary)
  } 
  # Case 3: order samples by their degree in the test network.
  else {
    # Order modules and samples by the test network
    sampleOrder <- sampleOrderInternal(testProps, verbose, TRUE)[[di]][[ti]][[1]]
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  
  # Case 1: di == ti. Plotting within the same dataset => nothing missing.
  # Case 2: orderBy == di: those missing in the discovery should have grey bars.
  # Case 3: orderBy == ti: ordering within the same dataset => nothing missing.
  
  if (orderNodesBy == "discovery" & di != ti) {
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

  if (orderSamplesBy == "discovery" & di != ti) {
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
  # Set up other property vectors and datasets
  #-----------------------------------------------------------------------------
  
  testProps <- testProps[[di]][[ti]] # collapse for easy access
  
  # (Normalised) weighted degree vector
  wDegreeVec <- foreach(mi = seq_along(testProps), .combine=c) %do% {
    testProps[[mi]]$degree/max(na.omit(testProps[[mi]]$degree))
  }
  wDegreeVec <- wDegreeVec[nodeOrder]
  
  if (!is.null(scaledData[[ti]])) {
    # node contribution
    nodeContribVec <- foreach(mi = seq_along(testProps), .combine=c) %do% {
      testProps[[mi]]$contribution
    }
    nodeContribVec <- nodeContribVec[nodeOrder]
    
    # and respective colors
    nodeContribCols <- rep(correlation.palette()[1], length(nodeContribVec))
    nodeContribCols[nodeContribVec > 0] <- tail(correlation.palette(), 1)
    
    # Summary profile matrix
    summaries <- foreach(mi = seq_along(testProps), .combine=cbind) %do% {
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
  }


  #-----------------------------------------------------------------------------
  # Set up plotting region
  #-----------------------------------------------------------------------------
  gaxt <- NULL
  if (plotNodeNames)
    gaxt <- nodeOrder
  
  if (missing(maxt.line) && !plotNodeNames) {
    maxt.line <- -0.5
  } else if (missing(maxt.line) && plotNodeNames) {
    maxt.line <- 3
  }
  
  if (is.null(scaledData[[ti]])) {
    # set up plot layout
    layout(
      mat=matrix(1:3, ncol=1), heights=c(0.4, 0.4, 0.2)
    )
  } else {
    # set up plot layout
    summary.window <- min(0.2 + (length(mods) - 1)* 0.1, 0.5) 
    layout(
      mat=matrix(c(rep(8, 5), 7, 1:6), ncol=2), 
      heights=c(0.7/3, 0.7/3, 0.12, 0.12, 0.06, 0.7/3),
      widths=c(summary.window, 1-summary.window)
    )
  }
  par(oma=old.par[["mar"]]+old.par[["oma"]])
  
  #-----------------------------------------------------------------------------
  # Plot topology components
  #-----------------------------------------------------------------------------
  # Plot correlation
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    correlation[[ti]][presentNodes, presentNodes], 
    correlation.palette(), c(-1, 1), 
    moduleAssignments[[di]][nodeOrder], na.pos.x, plotLegend=TRUE, 
    main="", legend.main="Correlation", plotModuleNames=FALSE, 
    legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
    legend.line=0.1, maxt.line=maxt.line, border.width=border.width
  )
  mtext(main, side=3, line=1, cex=par('cex.main'), font=2, xpd=NA)
  
  # Plot network
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    network[[ti]][presentNodes, presentNodes], network.palette(), 
    c(0, 1), moduleAssignments[[di]][nodeOrder], na.pos.x, 
    plotLegend=TRUE, main="", legend.main="Edge weights", 
    plotModuleNames=FALSE, legend.tick.size=legend.tick.size, 
    laxt.line=laxt.line, legend.line=0.1, maxt.line=maxt.line,
    border.width=border.width
  )
  
  # Plot weighted degree
  par(mar=c(1,1,1,1))
  if (is.null(scaledData[[ti]])) {
    plotBar(
      wDegreeVec, c(0,1), moduleAssignments[[di]][nodeOrder], "#feb24c", 
      drawBorders=drawBorders, plotModuleNames=plotModuleNames, 
      xaxt=plotNodeNames, xaxt.line=gaxt.line, main="",
      ylab="Weighted\ndegree", maxt.line=maxt.line, 
      border.width=border.width
    )
  } else {
    plotBar(
      wDegreeVec, c(0,1), moduleAssignments[[di]][nodeOrder], "#feb24c", 
      drawBorders=drawBorders, plotModuleNames=FALSE, main="", xaxt=FALSE,
      ylab="Weighted\ndegree", maxt.line=maxt.line, 
      border.width=border.width
    )
  }
  
  if (!is.null(scaledData[[ti]])) {
    # Plot Module Membership
    par(mar=c(1, 1, 1, 1))
    plotBar(
      nodeContribVec, c(-1,1), moduleAssignments[[di]][nodeOrder], 
      nodeContribCols, drawBorders=drawBorders, plotModuleNames=FALSE, main="", 
      xaxt=FALSE, ylab="Node\ncontribution", maxt.line=maxt.line, 
      border.width=border.width
    )
    
    # Plot the data matrix
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
    
    par(mar=c(1,1,1,1))
    emptyPlot(c(0,1), c(0,1), bty="n")

    yaxt <- NULL
    if (plotSampleNames)
      yaxt <- sampleOrder
    par(mar=c(1, 1, 1, 1))
    plotSquareHeatmap(
      dat, palette, vlim=range.pal, legend.lim=range.dat,
      moduleAssignments[[di]][nodeOrder], na.pos.x, na.pos.y, 
      xaxt=gaxt, yaxt=NULL, plotLegend=FALSE, main="",
      legend.main="", plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, maxt.line=maxt.line, border.width=border.width
    )
    
    # Plot data legend
    nNodes <- ncol(dat) + length(na.pos.x)
    nSamples <- nrow(dat) + length(na.pos.y)
    if (all(range.dat < 0)) {
      addGradientLegend(
        head(palette, length(palette)/2), range.dat, range.dat, TRUE, 
        xlim=c(0.5+nNodes*0.1,nNodes+0.5-nNodes*0.1), 
        ylim=c(nSamples+0.5+nSamples*0.2,nSamples+0.5+nSamples*0.3), 
        tick.size=legend.tick.size, border.width=border.width,
        main="Module data", axis.line=laxt.line, srt=0
      )
    } else if (all(range.dat > 0)) {
      addGradientLegend(
        tail(palette, length(palette)/2), range.dat, range.dat, TRUE, 
        xlim=c(0.5+nNodes*0.1,nNodes+0.5-nNodes*0.1), 
        ylim=c(nSamples+0.5+nSamples*0.2,nSamples+0.5+nSamples*0.3), 
        tick.size=legend.tick.size, border.width=border.width,
        main="Module data", axis.line=laxt.line, srt=0
      )
    } else {
      plim <- c(-max(abs(range.dat)), max(abs(range.dat)))
      addGradientLegend(
        palette, plim, range.dat, TRUE, main="Module data",
        xlim=c(0.5+nNodes*0.1,nNodes+0.5-nNodes*0.1), 
        ylim=c(nSamples+0.5+nSamples*0.2,nSamples+0.5+nSamples*0.3),  
        tick.size=legend.tick.size, border.width=border.width,
        axis.line=laxt.line, srt=0
      )
    }
    
    # Plot bar chart
    xlab <- "Module summary"
    if (length(mods) == 1) 
      xlab <- gsub(" ", "\n", xlab)
    par(mar=c(1, 1, 1, 1))
    plotMultiBar(
      summaries, rep(list(range(summaries, na.rm=TRUE)), ncol(summaries)),
      cols=summaries.cols , drawBorders=drawBorders, border.width=border.width,
      yaxt=plotSampleNames, plotModuleNames=plotModuleNames, 
      yaxt.line=saxt.line, maxt.line=0, xlab=xlab, 
      cex.modules=par("cex.lab")*0.7
    )
  }
}