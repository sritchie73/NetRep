#' Plot the topology of a network module
#' 
#' @description 
#' Plot the correlation structure, network edges, (normalised) connectivity, 
#' module membership, underlying data, and module summary vector of one or
#' more network modules.
#' 
#' Individual components of the module plot can be plotted using 
#' \code{\link{plotCorrelation}}, \code{\link{plotNetwork}}, 
#' \code{\link{plotConnectivity}}, \code{\link{plotModuleMembership}}, 
#' \code{\link{plotData}}, and \code{\link{plotModuleSummary}}.
#' 
#' @inheritParams common_params
#' @inheritParams common_params2
#' 
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
#' @param main title for the plot.
#' @param drawBorders logical; if \code{TRUE}, borders are drawn around the 
#'  connectivity, module membership, and module summary bar plots.
#' @param border.width line width for borders.
#' @param gaxt.line the number of lines into the bottom margin at which the node
#'  names will be drawn.
#' @param saxt.line the number of lines into the left margin at which the sample
#'  names will be drawn.
#' @param maxt.line the number of lines into the bottom margin at which the 
#'  module names will be drawn.
#' @param legend.tick.size size of the ticks on each axis legend relative to the
#'  size of the correlation, edge weights, and data matrix heatmaps.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{legend.tick.size}.
#' @param cex.axis relative size of the node and sample names.
#' @param cex.lab relative size of the module names, legend titles, and axis
#'  labels.
#' @param cex.main relative size of the plot title.
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
#' \code{\link{plotConnectivity}},
#' \code{\link{plotModuleMembership}},
#' \code{\link{plotData}}, and
#' \code{\link{plotModuleSummary}}.
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
  data, correlation, network, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderNodesBy="discovery",
  orderModules, plotNodeNames=TRUE, plotSampleNames=TRUE, plotModuleNames,
  main="Module Topology", drawBorders=FALSE, border.width=2, gaxt.line=-0.5, 
  saxt.line=-0.5, maxt.line, legend.tick.size=0.04, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
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
      if (!is.null(data[[discovery]])) {
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
      } else {
        warning(
          "Data matrix missing from the 'discovery' dataset, modules will ",
          "be ordered as provided"
        )
      }
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
      if (!is.null(data[[test]])) {
        seps <- matrix(
          0, ncol=length(props), 
          nrow=length(props[[1]]$moduleSummary)
        )
        colnames(seps) <- names(props)
        for (mi in seq_along(props)) {
          seps[,mi] <- props[[mi]]$moduleSummary
        }
        moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
      } else {
        warning(
          "Data matrix missing from the 'test' dataset, modules will be ",
          "ordered as provided"
        )
      }
    }
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  # If test == discovery, we order samples and modules by discovery, and plot
  # the discovery. So we need to make sure the data matrix exists in the
  # discovery dataset
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
  } else if (!is.null(data[[test]])) {
    if (orderSamplesBy == "none") {
      sampleOrder <- rownames(data[[test]])
    } else {
      sampleOrder <- names(sort(
        props[moduleOrder][[1]]$moduleSummary, decreasing=TRUE
      ))
    }
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  if (!is.null(data[[test]])) {
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
  }

  # Handle nodes not present in the test dataset
  na.pos.x <- which(nodeOrder %nin% colnames(correlation[[test]]))
  if (length(na.pos.x) > 0) {
    presentVars <- nodeOrder[-na.pos.x]
  } else {
    presentVars <- nodeOrder
  }
  
  #-----------------------------------------------------------------------------
  # Set up other property vectors and datasets
  #-----------------------------------------------------------------------------
  # (Normalised) Intramodular Connectivity vector
  kIM <- foreach(mi = seq_along(props), .combine=c) %do% {
    # Normalise the connectivity by the maximum. The value has no meaning,
    # just the relative sizes and ranks
    props[[mi]]$connectivity/max(na.omit(props[[mi]]$connectivity))
  }
  kIM <- kIM[nodeOrder]
  
  if (!is.null(data[[test]])) {
    # now build the Module Membership vector
    MM <- foreach(mi = seq_along(props), .combine=c) %do% {
      props[[mi]]$moduleMembership
    }
    MM <- MM[nodeOrder]
    MM.cols <- rep(correlation.palette()[1], length(MM))
    MM.cols[MM > 0] <- tail(correlation.palette(), 1)
    
    # Summary Expression profiles matrix
    SEP <- foreach(mi = moduleOrder, .combine=cbind) %do% {
      matrix(
        insert.nas(props[[mi]]$moduleSummary[presentSamples], na.pos.y),
        ncol=1
      )
    }
    colnames(SEP) <- moduleOrder
    rownames(SEP) <- sampleOrder
    
    # Now build the colors
    cols <- matrix(data.palette()[1], nrow(SEP), ncol(SEP))
    cols[SEP > 0] <- tail(data.palette(), 1)
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
  
  if (is.null(data[[test]])) {
    # set up plot layout
    layout(
      mat=matrix(1:3, ncol=1), heights=c(0.4, 0.4, 0.2)
    )
  } else {
    # set up plot layout
    SEPw <- min(0.2 + (length(modules) - 1)* 0.1, 0.5) 
    layout(
      mat=matrix(c(rep(8, 5), 7, 1:6), ncol=2), 
      heights=c(0.7/3, 0.7/3, 0.12, 0.12, 0.06, 0.7/3),
      widths=c(SEPw, 1-SEPw)
    )
  }
  par(oma=old.par[["mar"]]+old.par[["oma"]])
  
  #-----------------------------------------------------------------------------
  # Plot topology components
  #-----------------------------------------------------------------------------
  # Plot correlation
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    correlation[[test]][presentVars , presentVars], 
    correlation.palette(), c(-1, 1), 
    moduleAssignments[[discovery]][nodeOrder], na.pos.x, plotLegend=TRUE, 
    main="", legend.main="Correlation", plotModuleNames=FALSE, 
    legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
    legend.line=0.1, maxt.line=maxt.line, border.width=border.width
  )
  mtext(main, side=3, line=1, cex=par('cex.main'), font=2, xpd=NA)
  
  # Plot network
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    network[[test]][presentVars , presentVars], network.palette(), 
    c(0, 1), moduleAssignments[[discovery]][nodeOrder], na.pos.x, 
    plotLegend=TRUE, main="", legend.main="Edge weights", 
    plotModuleNames=FALSE, legend.tick.size=legend.tick.size, 
    laxt.line=laxt.line, legend.line=0.1, maxt.line=maxt.line,
    border.width=border.width
  )
  
  # Plot Intamodular Connectivity
  par(mar=c(1,1,1,1))
  if (is.null(data[[test]])) {
    plotBar(
      kIM, c(0,1), moduleAssignments[[discovery]][nodeOrder], "#feb24c", 
      drawBorders=drawBorders, plotModuleNames=plotModuleNames, 
      xaxt=plotNodeNames, xaxt.line=gaxt.line, main="",
      ylab="Normalised connectivity", maxt.line=maxt.line, 
      border.width=border.width
    )
  } else {
    plotBar(
      kIM, c(0,1), moduleAssignments[[discovery]][nodeOrder], "#feb24c", 
      drawBorders=drawBorders, plotModuleNames=FALSE, main="", xaxt=FALSE,
      ylab="Normalised\nconnectivity", maxt.line=maxt.line, 
      border.width=border.width
    )
  }
  
  if (!is.null(data[[test]])) {
    # Plot Module Membership
    par(mar=c(1, 1, 1, 1))
    plotBar(
      MM, c(-1,1), moduleAssignments[[discovery]][nodeOrder], MM.cols,
      drawBorders=drawBorders, plotModuleNames=FALSE, main="", xaxt=FALSE,
      ylab="Module\nMembership", maxt.line=maxt.line, border.width=border.width
    )
    
    # Plot the data matrix
    ge <- data[[test]][presentSamples, presentVars]
    range.dat <- range(ge)
    if (all(range.dat > 0)) {
      palette <- tail(data.palette(), length(data.palette())/2)
    } else if (all(range.dat < 0)) {
      palette <- head(data.palette(), length(data.palette())/2)
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
      ge, palette, vlim=range.pal, legend.lim=range.dat,
      moduleAssignments[[discovery]][nodeOrder], na.pos.x, na.pos.y, 
      xaxt=gaxt, yaxt=NULL, plotLegend=FALSE, main="",
      legend.main="", plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, maxt.line=maxt.line, border.width=border.width
    )
    
    # Plot data legend
    if (all(range.dat < 0)) {
      addGradientLegend(
        head(palette, length(palette)/2), range.dat, range.dat, TRUE, 
        xlim=c(0.5+ncol(ge)*0.1,ncol(ge)+0.5-ncol(ge)*0.1), 
        ylim=c(nrow(ge)+0.5+nrow(ge)*0.2,nrow(ge)+0.5+nrow(ge)*0.3), 
        tick.size=legend.tick.size, border.width=border.width,
        main="Module data", axis.line=laxt.line
      )
    } else if (all(range.dat > 0)) {
      addGradientLegend(
        tail(palette, length(palette)/2), range.dat, range.dat, TRUE, 
        xlim=c(0.5+ncol(ge)*0.1,ncol(ge)+0.5-ncol(ge)*0.1), 
        ylim=c(nrow(ge)+0.5+nrow(ge)*0.2,nrow(ge)+0.5+nrow(ge)*0.3), 
        tick.size=legend.tick.size, border.width=border.width,
        main="Module data", axis.line=laxt.line
      )
    } else {
      plim <- c(-max(abs(range.dat)), max(abs(range.dat)))
      addGradientLegend(
        palette, plim, range.dat, TRUE, main="Module data",
        xlim=c(0.5+ncol(ge)*0.1,ncol(ge)+0.5-ncol(ge)*0.1), 
        ylim=c(nrow(ge)+0.5+nrow(ge)*0.2,nrow(ge)+0.5+nrow(ge)*0.3),  
        tick.size=legend.tick.size, border.width=border.width,
        axis.line=laxt.line
      )
    }
    
    # Plot bar chart
    xlab <- "Module summary"
    if (length(modules) == 1) 
      xlab <- gsub(" ", "\n", xlab)
    par(mar=c(1, 1, 1, 1))
    plotMultiBar(
      SEP, rep(list(range(SEP, na.rm=TRUE)), ncol(SEP)),
      cols=cols, drawBorders=drawBorders, border.width=border.width,
      yaxt=plotSampleNames, plotModuleNames=plotModuleNames, 
      yaxt.line=saxt.line, maxt.line=0, xlab=xlab, 
      cex.modules=par("cex.lab")*0.7
    )
  }
}