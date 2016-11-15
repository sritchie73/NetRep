#' Map a variable to a color gradient
#'
#' For a color gradient mapped onto a scale of values (\code{vlim}), retrieve
#' the color for any given \code{value}
#'
#' @param values a vector of values to retrieve the colors for.
#' @param palette a vector of colors
#' @param vlim limits of the values
#'
#' @return
#'  The color for a given value.
#'  
#' @import utils
#' 
#' @keywords internal
getColFromPalette <- function(values, palette, vlim) {
  if (missing(vlim)) {
    vlim <- range(values)
  }
  vpal <- seq(vlim[1], vlim[2], length=length(palette)+1)

  sapply(values, function(vv) {
    if (vv >= vlim[2]) {
      tail(palette, n=1)
    } else if (vv <= vlim[1]) {
      palette[1]
    } else {
      palette[sum(vv > vpal)]
    }
  })
}

#' Get the plot limits to set for the desired plot window
#'
#' \code{plot} manually adds an extra region to the plot on top of the given
#' 'xlim' and 'ylim', amounting to 4% of the plot region in either direction.
#' This function tells you what limits to set so that the boundaries of the plot
#' are precisely the min and max you want.
#'
#' @param dlim the limits you want to set
#' @return
#'  the min and max to set for the desired plot limits
#'
#' @keywords internal
forceLim <- function(dlim) {
  A = matrix(c(1.04, -0.04, -0.04, 1.04), nrow=2, ncol=2)
  B = matrix(dlim, nrow=2, ncol=1)
  as.vector(solve(A, B))
}

#' Create an empty plot window with the specified limits
#'
#' @param ... extra arguments to pass to \code{\link[graphics]{plot}}.
#' @param xlim limits for the x axis.
#' @param ylim limits for the y axis.
#' @param xlab label for the x axis.
#' @param ylab label for the y axis.
#' @param hardlim logical; if \code{TRUE}, the plot window is exactly
#'  constrained to the specified \code{xlim} and \code{ylim}.
#'
#' @keywords internal
emptyPlot <- function(..., xlim, ylim, xlab="", ylab="", hardlim=TRUE) {
  if (!missing(xlim) && !missing(ylim)) {
    if(hardlim) {
      xlim <- forceLim(xlim)
      ylim <- forceLim(ylim)
    }
    plot(
      ..., 0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab=xlab,
      ylab=ylab
    )
  } else if (!missing(ylim)) {
    if (hardlim) {
      ylim <- forceLim(ylim)
    }
    plot(
      ..., ylim=ylim, type='n', xaxt='n', yaxt='n', xlab=xlab, ylab=ylab
    )
  } else if (!missing(xlim)) {
    if (hardlim) {
      xlim <- forceLim(xlim)
    }
    plot(
      ..., xlim=xlim, type='n', xaxt='n', yaxt='n', xlab=xlab, ylab=ylab
    )
  } else {
    plot(
      ..., type='n', xaxt='n', yaxt='n', xlab=xlab, ylab=ylab
    )
  }
}

#' Color palette for correlation heatmaps
#'
#' RColorBrewer palette "RdYlBu" with the middle color replaced with white.
#' This gives a nicer contrast than the "RdBu" palette
#' 
#' @import RColorBrewer
#' @keywords internal
correlation.palette <- function() {
  cols <- rev(brewer.pal(11, "RdYlBu"))
  cols[6] <- "#FFFFFF"
  cols
}

#' Color palette for network heatmaps
#'
#' RColorBrewer palette "RdYlBu" with the middle color replaced with white.
#' 
#' @keywords internal
network.palette <- function() {
  correlation.palette()[6:11]
}

#' Get module break points on the x-axis
#'
#' @param mas ordered subset of the moduleAssignments vector
#'
#' @return
#'  a vector of positions on the x-axis where one module begins and another ends
#'  
#' @keywords internal
getModuleBreaks <- function(mas) {
  sizes <- rle(mas)
  breaks = numeric(length(sizes$lengths) + 1)
  breaks[1] <- 0.5
  for (mi in seq_along(sizes$lengths)) {
    breaks[mi + 1] <- breaks[mi] + sizes$lengths[mi]
  }
  breaks
}

#' Get module mid-points on the x-axis
#'
#' @param mas ordered subset of the moduleAssignments vector
#'
#' @return
#'  a vector of positions on the x-axis indicating the centre of a module
#'  
#' @keywords internal
getModuleMidPoints <- function(mas) {
  breaks <- getModuleBreaks(mas)
  mids <- numeric(length(breaks) - 1)
  for (bi in seq_along(breaks)[-1]) {
    mids[bi - 1] <- (breaks[bi] - breaks[bi - 1])/2 + breaks[bi - 1]
  }
  mids
}

#' Check if a character vector of colors is valid
#' 
#' Courtesy of Josh O'Brien's stackoverflow answer at
#' \url{http://stackoverflow.com/a/13290832/2341679}
#' 
#' @param colvec a character vectors of colors (hex or name) to validate.
#' 
#' @keywords internal
areColors <- function(colvec) {
  sapply(colvec, function(col) {
    tryCatch(is.matrix(col2rgb(col)), error = function(e) FALSE)
  })
}

#' Get the network properties and order for a plot
#'
#' @param network list returned by \code{'processInput'}.
#' @param data data returned by \code{'processInput'}.
#' @param moduleAssignments list returned by \code{'processInput'}.
#' @param modules vector of modules to show on the plot.
#' @param di name of the discovery dataset.
#' @param ti name of the test dataset.
#' @param orderNodesBy vector returned by \code{'processInput'}.
#' @param orderSamplesBy vector returned by \code{'processInput'}.
#' @param orderModules vector returned by \code{'checkPlotArgs'}.
#' @param datasetNames vector returned by \code{'processInput'}.
#' @param nDatasets vector returned by \code{'processInput'}.
#' @param verbose logical; turn on verbose printing.
#' @param loadedIdx index of the currently loaded dataset.
#' @param dataLoaded currently loaded data matrix (may be NULL).
#' @param networkLoaded currently loaded network matrix.
#'
#' @keywords internal
plotProps <- function(
  network, data, moduleAssignments, modules, di, ti, orderNodesBy, 
  orderSamplesBy, orderModules, datasetNames, nDatasets, verbose, loadedIdx, 
  dataLoaded, networkLoaded
) {
  mods <- modules[[di]]
  mi <- NULL # suppresses CRAN note
  
  # Scenarios:
  # - No ordering of nodes + samples. We only need to calculate the network 
  #   properties for the 'test' dataset.
  # - Ordering of nodes only. We need to calculate the network properties in
  #   all datasets specified in 'orderNodesBy' (may be one or more) and in the
  #   'test' dataset (may or may not be specified in 'orderNodesBy').
  # - Ordering of samples only. We need to calculate the network properties in
  #   the 'orderSamplesBy' dataset, and in the 'test' dataset (which may or 
  #   may not be the same as 'orderSamplesBy').
  # - Ordering of both. We need to calculate the network properties in the
  #   'orderSamplesBy', 'orderNodesBy', and 'test' datasets.
  # this vector contains all datasets required for plotting
  
  plotDatasets <- unique(na.omit(c(orderSamplesBy, orderNodesBy, ti)))
  # The test dataset is the last plotDataset: this is so we don't have to
  # load it again after calculating the network properties
  plotDatasets <- c(plotDatasets[-which(plotDatasets == ti)], ti)
  plotDatasets <- list(plotDatasets)
  names(plotDatasets) <- datasetNames[di]
  
  # Calculate the network properties for all datasets required
  res <- netPropsInternal(network, data, moduleAssignments, modules, di, 
                          plotDatasets, nDatasets, datasetNames, verbose, 
                          loadedIdx, dataLoaded, networkLoaded, TRUE)
  props <- res$props
  loadedIdx <- res$loadedIdx
  
  # Order nodes based on degree
  if (length(orderNodesBy) > 1 || !is.na(orderNodesBy)) {
    if (length(orderNodesBy) > 1) {
      mean <- TRUE
    } else {
      mean <- FALSE
    }
    
    # nodeOrderInternal will average acros all test datasets, so we need to 
    # filter just to those specified in 'orderNodesBy' while preserving the
    # structure of 'props'
    orderProps <- filterInternalProps(props, orderNodesBy, di)
    nodeOrder <- nodeOrderInternal(
      orderProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE, mean
    )
    nodeOrder <- simplifyList(nodeOrder, depth=3)
    
    # The module order will be the names of the simplified list iff there are
    # multiple modules to render
    if (!is.list(nodeOrder)) {
      moduleOrder <- mods
      if (is.numeric(moduleOrder))
        moduleOrder <- as.character(moduleOrder)
    } else {
      moduleOrder <- names(nodeOrder)
    }
    
    # Now flatten the node order list
    nodeOrder <- unlist(nodeOrder)
  } else {
    hasProps <- !sapply(props[[di]][[ti]], is.null) 
    moduleOrder <- names(props[[di]][[ti]])[hasProps]
    nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[di]][[ti]][[mi]]$degree)
    }
  }
  
  if (is.null(orderSamplesBy)) {
    sampleOrder <- NULL
  } else if (!is.na(orderSamplesBy)) {
    orderProps <- filterInternalProps(props, orderSamplesBy, di, moduleOrder[1])
    sampleOrder <- sampleOrderInternal(orderProps, verbose, na.rm=FALSE)
    sampleOrder <- simplifyList(sampleOrder, depth=3)
  } else {
    sampleOrder <- rownames(deref(dataLoaded))
  }
  
  # Just keep the properties we need for plotting
  testProps <- simplifyList(props[[di]][[ti]], depth=1)
  if (length(moduleOrder) == 1) {
    testProps <- list(testProps)
    names(testProps) <- moduleOrder
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  
  na.pos.x <- which(nodeOrder %nin% colnames(deref(networkLoaded)))
  if (length(na.pos.x) > 0) {
    presentNodes <- nodeOrder[-na.pos.x]
  } else {
    presentNodes <- nodeOrder
  }
  
  if (is.null(sampleOrder)) {
    na.pos.y <- NULL
    presentSamples <- NULL
  } else if (!is.numeric(sampleOrder)) {
    na.pos.y <- which(sampleOrder %nin% rownames(deref(dataLoaded)))
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
  # Filter nodes and samples to display on the plot
  #-----------------------------------------------------------------------------
  # Unless we are ordering nodes/samples by the discovery dataset, we should 
  # only show nodes or samples present in the current dataset.
  if (length(orderNodesBy) == 1 && (is.na(orderNodesBy) || orderNodesBy != di)) {
    nodeOrder <- presentNodes
    na.pos.x <- NULL
  }
  
  if (!is.null(sampleOrder) && (is.na(orderSamplesBy) || orderSamplesBy != di)) {
    sampleOrder <- presentSamples
    na.pos.y <- NULL
  }
  
  
  # Returned processed results
  return(list(
    testProps=testProps, nodeOrder=nodeOrder, moduleOrder=moduleOrder,
    sampleOrder=sampleOrder, na.pos.x=na.pos.x, na.pos.y=na.pos.y,
    presentNodes=presentNodes, presentSamples=presentSamples
  ))
}
