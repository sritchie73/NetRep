#' Plot a symmetric heatmap as a triangle
#' 
#' @param values values to plot on the heatmap.
#' @param palette color palette to interpolate over.
#' @param vlim range of values to use when mapping values to the \code{palette}.
#' @param mas ordered subset of the moduleAssignments vector.
#' @param na.indices indices of missing values on the x axis.
#' @param na.col color of missing values to plot.
#' @param xaxt character vector of names to print along the x axis.
#' @param plotModuleNames logical; if \code{TRUE} the names of the modules are
#'  plotted along the x axis.
#' @param main title for the plot.
#' @param main.line the number of lines into the top margin at which the plot
#'  title will be drawn.
#' @param plotLegend logical; if \code{TRUE} a legend is added to the right 
#'  side of the plot.
#' @param legend.main title for the legend.
#' @param xaxt.line the number of lines into the margin at which the x axis 
#'  labels will be drawn.
#' @param maxt.line the number of lines into the margin at which the module 
#'  names will be drawn.
#' @param laxt.tck size of the ticks on the axis legend as a proportion
#'  of the horizontal size of the plot window.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{laxt.tck}.
#' @param legend.line the distance from the left of the plot to render the 
#'  legend as a proportion of the horizontal size of the plot window.
#' @param lwd line width for borders.
#' @param dryRun logical; if \code{TRUE} only the axes and borders will be 
#'  drawn.
#' 
plotTriangleHeatmap <- function(
  values, palette, vlim, mas, na.indices=NULL, na.col="#bdbdbd", xaxt=NULL,
  plotModuleNames=TRUE, main="", main.line=0, plotLegend=TRUE, legend.main="",
  xaxt.line=-0.5, maxt.line=3, laxt.tck=0.04, laxt.line=2.5, 
  legend.line=0.1, lwd=2, dryRun=FALSE
) {
  nNodes <- ncol(values) + length(na.indices)
  palette <- colorRampPalette(palette)(255)
  
  if (length(vlim) == 1) {
    vlim <- c(0, max(values[lower.tri(values)]))
  }
  
  # Use a fixed width/height for all plots so that offsets and margin lines
  # are the same regardless of the number of nodes shown
  unitSize <- 1/nNodes
  halfUnit <- unitSize/2
  ph <- (nNodes * unitSize)/2
  pw <- nNodes * unitSize
  
  # Create empty plot
  emptyPlot(xlim=c(halfUnit, pw + halfUnit), ylim=c(0, ph), bty="n")
  
  # render triangles row by row
  if (!dryRun) {
    for (plotRow in 1:nNodes) {
      startCol <- nNodes - (plotRow - 1)
      for (ii in 1:plotRow) {
        jj <- startCol + (ii - 1)
        ci <- ii
        cj <- jj
        
        topy <- (nNodes - (plotRow - 1))/2
        # If we're on the diagonal, plot a triangle, otherwise a diamond
        if (plotRow == nNodes) {
          boty <- 0 
        } else {
          boty <- topy - 1
        }
        
        xOffset <- (nNodes - (plotRow - 1))/2 
        rightx <- ii + xOffset
        leftx <- rightx - 1
        
        if (ii %nin% na.indices && jj %nin% na.indices) {
          col <- getColFromPalette(values[ci, cj], palette, vlim)
          cj <- cj + 1
          ci <- ci + 1
        } else {
          col <- na.col
        }
        
        leftx <- leftx * unitSize
        rightx <- rightx * unitSize
        boty <- boty * unitSize
        topy <- topy * unitSize
        
        polygon(
          x=c(leftx, leftx+halfUnit, rightx, leftx+halfUnit, leftx),
          y=c(topy-halfUnit, topy, topy-halfUnit, boty, topy-halfUnit),
          col=col, border=col
        )
      }
    } 
  }
  
  # render module boundaries
  if (length(unique(mas)) > 1) {
    breaks <- getModuleBreaks(mas)
    mids <- getModuleMidPoints(mas)
    for (mi in seq_along(mids)) {
      height <- (breaks[mi + 1] - mids[mi]) * unitSize
      leftx <- breaks[mi] * unitSize
      rightx <- breaks[mi+1] * unitSize
      midx <- mids[mi] * unitSize
      polygon(
        x=c(leftx, rightx, midx, leftx),
        y=c(0, 0, height, 0), lwd=lwd
      )
    }
  }
  if (plotModuleNames) {
    axis(
      side=1, las=1, 
      at=getModuleMidPoints(mas)*unitSize,
      labels=unique(mas), line=maxt.line, tick=FALSE,
      cex.axis=par("cex.lab"), font=2
    )
  }
  
  # render border of plot
  polygon(
    x=c(halfUnit, pw + halfUnit, ph + halfUnit, halfUnit), 
    y=c(0, 0, ph, 0),
    lwd=lwd, xpd=NA
  )
  
  # Render axes
  if (!is.null(xaxt)) {
    axis(
      side=1, las=2, tick=FALSE, line=xaxt.line,
      at=(1:nNodes) * unitSize, labels=xaxt
    )
  }
  mtext(main, side=3, cex=par("cex")*par("cex.main"), font=2, line=main.line)
  
  # Add legend if specified
  if (plotLegend) {
    addGradientLegend(
      palette, vlim, TRUE, legend.main,
      xlim=c(halfUnit - pw*legend.line, pw*0.25), 
      ylim=c(ph/2 + ph*0.17, ph/2 + ph*0.25), tck=laxt.tck,
      axis.line=laxt.line, lwd=lwd
    )
  }
}

#' Plot a square heatmap
#' 
#' @param values values to plot on the heatmap.
#' @param palette color palette to interpolate over.
#' @param vlim range of values to use when mapping values to the \code{palette}.
#' @param mas ordered subset of the moduleAssignments vector.
#' @param na.indices.x indices of missing values on the x axis.
#' @param na.indices.y indices of missing values on the y axis.
#' @param na.col color of missing values to plot.
#' @param xaxt character vector of names to print along the x axis.
#' @param yaxt character vector of names to print along the y axis.
#' @param plotModuleNames logical; if \code{TRUE} the names of the modules are
#'  plotted along the x axis if \code{values} is not symmetric, and along both
#'  axes if \code{values} is symettric.
#' @param main title for the plot.
#' @param main.line the number of lines into the top margin at which the plot
#'  title will be drawn.
#' @param plotLegend logical; if \code{TRUE} a legend is added to the right 
#'  side of the plot.
#' @param legend.main title for the legend.
#' @param xaxt.line the number of lines into the margin at which the x axis 
#'  labels will be drawn.
#' @param yaxt.line the number of lines into the margin at which the y axis 
#'  labels will be drawn.
#' @param maxt.line the number of lines into the margin at which the module 
#'  names will be drawn.
#' @param laxt.tck size of the ticks on the axis legend as a proportion
#'  of the horizontal size of the plot window.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{laxt.tck}.
#' @param legend.line the distance from the plot to render the legend as a 
#'  proportion of the horizontal size of the plot window.
#' @param lwd line width for borders.
#' @param dryRun logical; if \code{TRUE} only the axes and borders will be 
#'  drawn.
#'  
plotSquareHeatmap <- function(
  values, palette, vlim, mas, na.indices.x=NULL, na.indices.y=NULL,
  na.col="#bdbdbd", xaxt=NULL, yaxt=NULL, plotModuleNames=TRUE, 
  main="", main.line=0, plotLegend=TRUE, legend.main="", xaxt.line=-0.5, 
  yaxt.line=-0.5, maxt.line=3, laxt.tck=0.04, laxt.line=2.5, 
  legend.line=0.1, lwd=2, dryRun=FALSE
) {
  nX <- ncol(values) + length(na.indices.x)
  nY <- nrow(values) + length(na.indices.y)
  palette <- colorRampPalette(palette)(255)
  
  if (length(vlim) < 2) {
    vlim <- c(0, max(c(values[lower.tri(values)], values[upper.tri(values)])))
  }
    
  # Use a fixed width/height for all plots so that offsets and margin lines
  # are the same regardless of the number of nodes shown
  xUnitSize <- 1/nX
  xHalfUnit <- xUnitSize/2
  yUnitSize <- 1/nY
  yHalfUnit <- yUnitSize/2
  
  pw <- nX * xUnitSize
  ph <- nY * yUnitSize
  
  # Create empty plot
  emptyPlot(xlim=c(xHalfUnit, pw + xHalfUnit), 
            ylim=c(yHalfUnit, ph + yHalfUnit), 
            bty="n")
  
  # render squares
  if (!dryRun) {
    cj <- 1
    for (jj in 1:nX) {
      ci <- 1
      for (ii in 1:nY) {
        if (ii %nin% na.indices.y && jj %nin% na.indices.x) {
          col <- getColFromPalette(values[ci, cj], palette, vlim)
          ci <- ci + 1
        } else {
          col <- na.col
        }
        
        xleft <- jj * xUnitSize - xHalfUnit
        xright <- jj * xUnitSize + xHalfUnit
        ybottom <- (nY - (ii - 1)) * yUnitSize - yHalfUnit
        ytop <- (nY - (ii - 1)) * yUnitSize + yHalfUnit
        
        rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col, 
             border=col)
      }
      if (jj %nin% na.indices.x) {
        cj <- cj + 1
      }
    }    
  }
  
  # render module boundaries
  if (length(unique(mas)) > 1) {
    breaks <- getModuleBreaks(mas)
    for (mi in seq_along(breaks)[-1]) {
      xleft <- breaks[mi - 1] * xUnitSize
      xright <- breaks[mi] * xUnitSize
      
      if (nX != nY) {
        ybottom <- yHalfUnit
        ytop <-  ph + yHalfUnit
      } else {
        ybottom <- (pw + xHalfUnit) - (breaks[mi] * yUnitSize - yHalfUnit)
        ytop <- (pw + xHalfUnit) - (breaks[mi - 1] * yUnitSize - yHalfUnit)
      }
      rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, 
           border="black", lwd=lwd)
    }
  }
  if (plotModuleNames) {
    if(!(nX == nY && is.null(xaxt) && !is.null(yaxt))) {
      axis(
        side=1, las=1, 
        at=getModuleMidPoints(mas) * xUnitSize,
        labels=unique(mas), line=maxt.line, tick=FALSE,
        cex.axis=par("cex.lab"), font=2
      )
    }
    if (nX == nY) {
      axis(
        side=2, las=2,
        at=ph + yHalfUnit - getModuleMidPoints(mas) * yUnitSize,
        labels=unique(mas), line=maxt.line, tick=FALSE,
        cex.axis=par("cex.lab"), font=2
      )
    }
  }
  
  # render border of plot
  rect(
    xleft=par("usr")[1],
    xright=par("usr")[2],
    ybottom=par("usr")[3],
    ytop=par("usr")[4],
    border="black",
    xpd=NA,
    lwd=lwd
  )
  
  # Render axes
  if (!is.null(xaxt)) {
    axis(
      side=1, las=2, tick=FALSE, line=xaxt.line,
      at=(1:nX)*xUnitSize, labels=xaxt
    )
  }
  if (!is.null(yaxt)) {
    axis(
      side=2, las=2, tick=FALSE, line=yaxt.line,
      at=(nY:1)*yUnitSize, labels=yaxt
    )
  }
  mtext(main, side=3, cex=par("cex")*par("cex.main"), font=2, line=main.line)
  
  # Add legend if specified
  if (plotLegend) {
    addGradientLegend(
      palette, vlim, FALSE, legend.main,
      xlim=c(pw - xHalfUnit + pw*legend.line, pw - xHalfUnit + pw*(legend.line+0.05)), 
      ylim=c(ph/3, ph - yHalfUnit - ph*0.1), tck=laxt.tck,
      axis.line=laxt.line, lwd=lwd
    )
  }
}

#' Plot a color palette legend
#' 
#' Add a legend to a plot window.
#' 
#' @param palette color palette.
#' @param legend.vlim limits of the values to display on the legend
#' @param horizontal logical; if \code{TRUE} the legend is plotted horizontally,
#'   otherwise vertically.
#' @param main title of the legend.
#' @param xlim xlim relative to the plotting region of the rest of the plot.
#' @param ylim ylim relative to the plotting region of the rest of the plot.
#' @param tck size of the legend axis ticks relative to the size of the 
#'  plot window.
#' @param axis.line distance from the axis to render the axis labels as a 
#'  multiple of \code{tck}.
#' @param lwd line width for borders.
#' @param srt angle of text labels
#' 
addGradientLegend <- function(
  palette, legend.vlim, horizontal, main, xlim, ylim, 
  tck=0.04, axis.line=3, lwd=2, srt
) {
  palette <- colorRampPalette(palette)(255)
  
  if (missing(xlim))
    xlim <- par("usr")[1:2]
  if (missing(ylim))
    ylim <- par("usr")[3:4]

  # Handle legends where the range of values doesn't map to the palette range 
  # (i.e. drawing a legend for the adjacency plot using the correlation 
  # palette)
  plim <- c(
    head(which(palette == getColFromPalette(
      legend.vlim[1], palette, legend.vlim
    )), 1),
    tail(which(palette == getColFromPalette(
      legend.vlim[2], palette, legend.vlim
    )), 1)
  )
  palette <- palette[plim[1]:plim[2]]
  # Draw gradient bar
  if (horizontal) {
    breaks <- seq(xlim[1], xlim[2], length=length(palette) + 1)
    for (pi in seq_along(palette)) {
      rect(
        xleft=breaks[pi],
        xright=breaks[pi + 1],
        ybottom=ylim[1],
        ytop=ylim[2],
        col=palette[pi],
        border=palette[pi],
        xpd=NA
      )
    }
  } else {
    breaks <- seq(ylim[1], ylim[2], length=length(palette) + 1)
    for (pi in seq_along(palette)) {
      rect(
        xleft=xlim[1],
        xright=xlim[2],
        ybottom=breaks[pi],
        ytop=breaks[pi+1],
        col=palette[pi],
        border=palette[pi],
        xpd=NA
      )
    }
  }
  
  # Render bounding box
  rect(
    xleft=xlim[1], xright=xlim[2], ybottom=ylim[1], ytop=ylim[2],
    border="black", lwd=lwd, xpd=NA
  )
  
  # Make sure axis ticks are centred at 0 if within the range of legend.vlim
  labels <- seq.int(legend.vlim[1L], legend.vlim[2L], length.out=5)
  labels <- format(labels, digits=2)
  if (horizontal) {
    tck <- (par("usr")[4] - par("usr")[3])*tck
    
    # for mapping from vlim to plot space
    v.per.x <- (xlim[2] - xlim[1])/(legend.vlim[2] - legend.vlim[1])
    zero <- xlim[1] +  v.per.x * (0 - legend.vlim[1])
    
    at <- seq.int(xlim[1L], xlim[2L], length.out=5)
    
    # Now plot the lines and text
    sapply(at, function(aa) {
      lines(x=c(aa, aa), y=c(ylim[1], ylim[1]-tck), lwd=lwd, xpd=NA)
    })
    
    if (missing(srt))
      srt <- 45
    
    text(labels, x=at, y=ylim[1]-tck*axis.line, cex=par("cex.axis"), xpd=NA, 
         srt=srt, adj=1)
  } else {
    tck <- (par("usr")[2] - par("usr")[1])*tck
    
    # for mapping from vlim to plot space
    v.per.y <- (ylim[2] - ylim[1])/(legend.vlim[2] - legend.vlim[1])
    zero <- ylim[1] +  v.per.y * (0 - legend.vlim[1])
    at <- seq.int(ylim[1L], ylim[2L], length.out=5)
    
    # draw axis ticks
    sapply(at, function(aa) {
      lines(x=c(xlim[1], xlim[1]-tck), y=c(aa, aa), lwd=lwd, xpd=NA)
    })
    if (missing(srt))
      srt <- 0
    text(labels, x=xlim[1]-tck*axis.line, y=at, cex=par("cex.axis"), xpd=NA)
  }
  
  # Render title
  offset <- (par("usr")[4] - par("usr")[3]) * 0.07
  text(
    main, x=xlim[1]+(xlim[2]-xlim[1])/2, y=ylim[2]+offset, font=2, xpd=NA,
    cex=par("cex.lab")
  )
}

#' Custom bar plot function
#' 
#' Plot bars around 0
#' 
#' @param heights heights of the bars.
#' @param heights.lim limits of the height axis.
#' @param mas ordered subset of the moduleAssignments vector
#' @param cols colors of each bar.
#' @param bar.width value between 0 and 1 controlling the proportion of space
#'  taken by each bar.
#' @param drawBorders logical; if \code{TRUE} a border is drawn around each bar.
#' @param na.col color of missing values to plot.
#' @param xaxt logical; If \code{TRUE}, the names of \code{heights} will be 
#'  rendered underneath the bar chart
#' @param plotModuleNames logical; if \code{TRUE} the names of the modules are
#'  plotted along the x axis.
#' @param main title for the plot.
#' @param main.line the number of lines into the top margin at which the plot
#'  title will be drawn.
#' @param xaxt.line the number of lines into the margin at which the x axis 
#'  labels will be drawn.
#' @param yaxt.line the number of lines into the margin at which the y axis
#'  tick labels will be drawn.
#' @param yaxt.tck the size of the y-axis tick marks.
#' @param maxt.line the number of lines into the margin at which the module 
#'  names will be drawn.
#' @param ylab label for the y axis
#' @param ylab.line the number of lines into the left margin at which the 
#'  y axis label will be drawn. 
#' @param lwd line width for borders.
#' @param dryRun logical; if \code{TRUE} only the axes and borders will be 
#'  drawn.
#'  
plotBar <- function(
  heights, heights.lim, mas, cols, bar.width=1, drawBorders=FALSE, 
  na.col="#bdbdbd", xaxt=TRUE, plotModuleNames=TRUE, main="", main.line=0,
  xaxt.line=-0.5, yaxt.line=0, yaxt.tck=-0.15, maxt.line=3, ylab="", 
  ylab.line=2.5, lwd=2, dryRun=FALSE
) {
  # Create vector of colors, one for each bar
  if (length(cols) == 1) {
    colvec <- rep(cols, length(heights))
  } else if (length(cols) == 2) { # Assume positive and negative
    colvec <- character(length=length(heights))
    colvec[heights > 0] <- cols[1]
    colvec[heights <= 0] <- cols[2]
  } else {
    stop("invalid length for 'cols'")
  }

  ylim <- heights.lim
  ylim[2] <- ylim[2] + (ylim[2] - ylim[1])*0.01
  ylim[1] <- ylim[1] - (ylim[2] - ylim[1])*0.01
  emptyPlot(
    xlim=c(0.5, length(heights)+0.5), ylim=ylim, bty="n", ylab="", xpd=NA
  )
  
  # Draw y-axis label
  mtext(ylab, side=2, cex=par("cex")*par("cex.lab"), font=1, line=ylab.line)
  
  # draw NAs
  if (!dryRun) {
    for (ii in seq_along(heights)) {
      if (is.na(heights[ii])) {
        rect(
          xleft=ii - 0.5,
          xright=ii + 0.5,
          ybottom=heights.lim[1],
          ytop=heights.lim[2],
          col=na.col,
          border=NA
        ) 
      }
    }
    
    for (ii in seq_along(heights)) {
      rect(
        xleft=ii-bar.width/2,
        xright=ii+bar.width/2,
        ybottom=0,
        ytop=heights[ii],
        col=colvec[ii],
        border=ifelse(drawBorders, "black", NA),
        lwd=lwd
      ) 
    }
  }
  abline(h=0, col="black", lwd=lwd)
  
  
  # draw axis
  axis(side=2, labels=FALSE, tck=yaxt.tck, lwd=lwd, at=axTicks(side=2))
  axis(side=2, tick=FALSE, line=yaxt.line, las=2, at=axTicks(side=2), 
       labels=axTicks(side=2))
  
  # render module boundaries
  if (length(unique(mas)) > 1) {
    breaks <- getModuleBreaks(mas)
    for (bi in head(breaks[-1], -1)) {
      lines(x=rep(bi, 2), y=heights.lim, lwd=lwd)
    }
  }
  
  if (plotModuleNames) {
    axis(
      side=1, las=1, 
      at=getModuleMidPoints(mas),
      labels=unique(mas), line=maxt.line, tick=FALSE,
      cex.axis=par("cex.lab"), font=2, xpd=NA
    )
  }
  
  # Render axes
  if (xaxt) {
    axis(
      side=1, las=2, tick=FALSE, line=xaxt.line,
      at=1:length(heights), labels=names(heights)
    )
  }
  mtext(main, side=3, cex=par("cex")*par("cex.main"), font=2, line=main.line)
}

#' Plot multiple horizontal bar plots
#' 
#' @param lengths a matrix whose columns contain the lengths of each bar for 
#'  the given property (e.g. each column should be a module, or phenotype).
#' @param lengths.lim a list of limits for the lengths axes.
#' @param cols a matrix of colors for each bar.
#' @param bar.width value between 0 and 1 controlling the proportion of space
#'  taken by each bar.
#' @param drawBorders logical; if \code{TRUE} a border is drawn around each bar.
#' @param main title for the plot
#' @param main.line the number of lines into the top margin at which the plot
#'  title will be drawn.  
#' @param na.col color of missing values to plot.
#' @param yaxt logical; If \code{TRUE}, the rownames of \code{heights} will be 
#'  rendered to the left of the bars.
#' @param plotModuleNames logical; if \code{TRUE} the names of the modules are
#'  plotted along the x axis.
#' @param yaxt.line the number of lines into the margin at which the y axis 
#'  labels will be drawn.
#' @param maxt.line the number of lines into the margin at which the module 
#'  labels will be drawn.
#' @param xaxt.line the number of lines into the margin at which the x axis
#'  labels will be drawn
#' @param xaxt.tck the size of the x-axis ticks.
#' @param xlab x axis label
#' @param xlab.line the number of lines into the bottom margin at which the 
#'  x axis label will be drawn. 
#' @param cex.modules relative size of module names.
#' @param lwd line width for borders.
#' @param dryRun logical; if \code{TRUE} only the axes and borders will be 
#'  drawn.
#'  
plotMultiBar <- function(
  lengths, lengths.lim, cols, bar.width=1, drawBorders=FALSE, main="", 
  main.line=1, na.col="#bdbdbd", yaxt=TRUE, plotModuleNames=TRUE, yaxt.line=0, 
  maxt.line=2.5, xaxt.line=0, xaxt.tck=-0.025, xlab="", xlab.line=2.5, 
  cex.modules=0.7, lwd=2, dryRun=FALSE
) {
  if (!is.matrix(lengths))
    lengths <- matrix(lengths, ncol=lengths)
  if (missing(lengths.lim)) {
    lengths.lim <- lapply(seq_len(ncol(lengths)), function(ci) {
      range(lengths[,ci])
    })
  }
  if (!is.list(lengths.lim))
    lengths.lim <- list(lengths.lim)
  
  # Create vector of colors, one for each bar
  if (length(cols) == 1) {
    colmat <- matrix(cols, nrow=nrow(lengths), ncol=ncol(lengths))
  } else if (length(cols) == 2) { # Assume positive and negative
    colmat <- matrix("", nrow=nrow(lengths), ncol=ncol(lengths))
    colmat[lengths > 0] <- cols[1]
    colmat[lengths <= 0] <- cols[2]
  } else {
    stop("invalid length for 'cols'")
  }
  
  pw <- 0.7 # width of each plot within the 0-1 space
  
  emptyPlot(
    xlim=c(0, ncol(lengths)), ylim=c(0, nrow(lengths)*1.01), bty="n", xlab="",
    xpd=NA
  )
  
  # Draw y-axis label
  mtext(xlab, side=1, cex=par("cex")*par("cex.lab"), font=1, line=xlab.line)
  
  # draw NAs
  if (!dryRun) {
    for (jj in seq_len(nrow(lengths))) {
      if (all(is.na(lengths[jj,]))) {
        rect(
          xleft=0,
          xright=ncol(lengths),
          ybottom=nrow(lengths) - jj,
          ytop=nrow(lengths) - (jj - 1),
          col=na.col,
          border=NA
        ) 
      }
    }
  }
  
  for (ii in seq_len(ncol(lengths))) {
    # we need to map from the value range to a range of 0-1
    rr <- lengths.lim[[ii]]
    rr.size <- rr[2] - rr[1]
    
    if (min(rr) > 0) {
      ax <- min(rr)
    } else if (max(rr) < 0) {
      ax <- max(rr)
    } else {
      ax <- 0
    }
    
    # Get x position for an value in range rr
    getX <- function(val) {
      (ii - 1) + (1-pw)/2 + pw/rr.size * (val - rr[1])
    }
    
    # Only draw bars if dryRun is FALSE
    if (!dryRun) {
      for (jj in seq_len(nrow(lengths))) {
        rect(
          xleft=getX(ax),
          xright=getX(lengths[jj,ii]),
          ybottom=nrow(lengths) - jj + (1 - bar.width)/2,
          ytop=nrow(lengths) - (jj - 1) - (1 - bar.width)/2,
          col=colmat[jj, ii],
          border=ifelse(drawBorders, "black", NA),
          lwd=lwd
        ) 
      }
    }
    
    # draw 0 axis
    lines(x=rep(getX(ax), 2), y=c(0, nrow(lengths)), lwd=lwd)
    
    # draw axis
    axis(
      side=1, labels=FALSE, tck=xaxt.tck, lwd=lwd,
      at=unique(c(getX(rr[1]), getX(ax), getX(rr[2])))
    )
    axis(
      side=1, tick=FALSE, line=xaxt.line, las=2,
      at=unique(c(getX(rr[1]), getX(ax), getX(rr[2]))), 
      labels=prettyNum(unique(c(rr[1], ax, rr[2])), digits=2)
    )
    if (plotModuleNames) {
      mtext(
        colnames(lengths)[ii], side=3, at=ii-0.5, cex=par("cex")*cex.modules, 
        font=2, line=maxt.line
      )
    }
  }
    
  # Draw title
  mtext(
    main, side=3, at=ncol(lengths)/2, cex=par("cex")*par("cex.main"), font=2, 
    adj=0.5, line=main.line
  )
   
  # Draw sample names
  if (yaxt) {
    axis(
      side=2, tick=FALSE, las=2, at=1:nrow(lengths)-0.5,
      labels=rownames(lengths), line=yaxt.line
    )
  }
}
