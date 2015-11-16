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
#' @param plotLegend logical; if \code{TRUE} a legend is added to the right side of
#'  the plot.
#' @param legend.main title for the legend.
#' @param legend.lim range of values to show on the legend.
#' @param xaxt.line the number of lines into the margin at which the x axis 
#'  labels will be drawn.
#' @param maxt.line the number of lines into the margin at which the module 
#'  names will be drawn.
#' @param legend.tick.size size of the ticks on the axis legend as a proportion
#'  of the horizontal size of the plot window.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{legend.tick.size}.
#' @param legend.line the distance from the left of the plot to render the 
#'  legend as a proportion of the horizontal size of the plot window.
#' @param border.width line width for borders.
#' 
plotTriangleHeatmap <- function(
  values, palette, vlim, mas, na.indices=NULL, na.col="#bdbdbd", xaxt=NULL,
  plotModuleNames=TRUE, main="", plotLegend=TRUE, legend.main="", legend.lim, 
  xaxt.line=-0.5, maxt.line=3, legend.tick.size=0.04, laxt.line=2.5, 
  legend.line=0.1, border.width=2
) {
  nGenes <- ncol(values) + length(na.indices)
  emptyPlot(xlim=c(0.5, nGenes + 0.5), ylim=c(0, nGenes/2), bty="n")
  palette <- colorRampPalette(palette)(255)
  
  # render squares / triangles
  ci <- 1
  for (ii in 1:nGenes) {
    cj <- 1
    for (jj in 1:ii) {
      plotRow <- ii - jj
      topy <- (plotRow + 1)/2
      # If we're on the diagonal, plot a triangle, otherwise a diamond
      if (plotRow == 0) {
        boty <- 0
      } else {
        boty <- topy - 1
      }
      xOffset <- plotRow/2
      rightx <- jj + xOffset + 0.5
      leftx <- rightx - 1
      
      if (ii %nin% na.indices && jj %nin% na.indices) {
        col <- getColFromPalette(values[ci, cj], palette, vlim)
        cj <- cj + 1
      } else {
        col <- na.col
      }
      polygon(
        x=c(leftx, leftx+0.5, rightx, leftx+0.5, leftx),
        y=c(topy-0.5, topy, topy-0.5, boty, topy-0.5),
        col=col, border=col
      )
    }
    if (ii %nin% na.indices) {
      ci <- ci + 1
    }
  }
  
  # render module boundaries
  if (length(unique(mas)) > 1) {
    breaks <- getModuleBreaks(mas)
    mids <- getModuleMidPoints(mas)
    for (mi in seq_along(mids)) {
      height <- breaks[mi + 1] - mids[mi]
      polygon(
        x=c(breaks[mi], breaks[mi+1], mids[mi], breaks[mi]),
        y=c(0, 0, height, 0), lwd=border.width
      )
    }
  }
  if (plotModuleNames) {
    axis(
      side=1, las=1, 
      at=getModuleMidPoints(mas),
      labels=unique(mas), line=maxt.line, tick=FALSE,
      cex.axis=par("cex.lab"), font=2
    )
  }
  
  # render border of plot
  polygon(
    x=c(0.5, nGenes+0.5, nGenes/2+0.5, 0.5), y=c(0, 0, nGenes/2, 0),
    lwd=border.width, xpd=NA
  )
  
  # Render axes
  if (!is.null(xaxt)) {
    axis(
      side=1, las=2, tick=FALSE, line=xaxt.line,
      at=1:nGenes, labels=xaxt
    )
  }
  mtext(main, side=3, cex=par("cex.main"), font=2)
  
  # Add legend if specified
  if (plotLegend) {
    if (missing(legend.lim))
      legend.lim <- vlim
    ph <- nGenes/2
    pw <- nGenes
    addGradientLegend(
      palette, vlim, legend.lim, TRUE, legend.main,
      xlim=c(0.5 - pw*legend.line, pw*0.27), 
      ylim=c(ph/2 + ph*0.17, ph/2 + ph*0.27), tick.size=legend.tick.size,
      axis.line=laxt.line, border.width=border.width
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
#' @param plotLegend logical; if \code{TRUE} a legend is added to the right side of
#'  the plot.
#' @param legend.main title for the legend.
#' @param legend.lim range of values to show on the legend.
#' @param xaxt.line the number of lines into the margin at which the x axis 
#'  labels will be drawn.
#' @param yaxt.line the number of lines into the margin at which the y axis 
#'  labels will be drawn.
#' @param maxt.line the number of lines into the margin at which the module 
#'  names will be drawn.
#' @param legend.tick.size size of the ticks on the axis legend as a proportion
#'  of the horizontal size of the plot window.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{legend.tick.size}.
#' @param legend.line the distance from the plot to render the legend as a 
#'  proportion of the horizontal size of the plot window.
#' @param border.width line width for borders.
#' 
plotSquareHeatmap <- function(
  values, palette, vlim, mas, na.indices.x=NULL, na.indices.y=NULL,
  na.col="#bdbdbd", xaxt=NULL, yaxt=NULL, plotModuleNames=TRUE, 
  main="", plotLegend=TRUE, legend.main="", legend.lim, xaxt.line=-0.5, 
  yaxt.line=-0.5, maxt.line=3, legend.tick.size=0.04, laxt.line=2.5, 
  legend.line=0.1, border.width=2
) {
  nX <- ncol(values) + length(na.indices.x)
  nY <- nrow(values) + length(na.indices.y)
  emptyPlot(xlim=c(0.5, nX+0.5), ylim=c(0.5, nY+0.5), bty="n")
  palette <- colorRampPalette(palette)(255)
  
  # render squares / triangles
  ci <- 1
  for (ii in 1:nY) {
    cj <- 1
    for (jj in 1:nX) {
      if (ii %nin% na.indices.y && jj %nin% na.indices.x) {
        col <- getColFromPalette(values[ci, cj], palette, vlim)
        cj <- cj + 1
      } else {
        col <- na.col
      }
      rect(
        xleft = jj - 0.5,
        xright = jj + 0.5,
        ybottom = (nY - (ii - 1)) - 0.5,
        ytop = (nY- (ii - 1)) + 0.5,
        col=col, border=col
      )
    }
    if (ii %nin% na.indices.y) {
      ci <- ci + 1
    }
  }
  
  # render module boundaries
  if (length(unique(mas)) > 1) {
    breaks <- getModuleBreaks(mas)
    for (mi in seq_along(breaks)[-1]) {
      if (nX != nY) {
        rect(
          xleft = breaks[mi - 1],
          xright = breaks[mi],
          ybottom = 0.5,
          ytop = nY + 0.5,
          border="black", lwd=border.width
        )
      } else {
        rect(
          xleft = breaks[mi - 1],
          xright = breaks[mi],
          ybottom = (nX + 0.5) - (breaks[mi] - 0.5),
          ytop = (nX + 0.5) - (breaks[mi - 1] - 0.5),
          border="black", lwd=border.width
        )
      }
    }
  }
  if (plotModuleNames) {
    if(!(nX == nY && is.null(xaxt) && !is.null(yaxt))) {
      axis(
        side=1, las=1, 
        at=getModuleMidPoints(mas),
        labels=unique(mas), line=maxt.line, tick=FALSE,
        cex.axis=par("cex.lab"), font=2
      )
    }
    if (nX == nY) {
      axis(
        side=2, las=2,
        at=nY + 0.5 - getModuleMidPoints(mas),
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
    lwd=border.width
  )
  
  # Render axes
  if (!is.null(xaxt)) {
    axis(
      side=1, las=2, tick=FALSE, line=xaxt.line,
      at=1:nX, labels=xaxt
    )
  }
  if (!is.null(yaxt)) {
    axis(
      side=2, las=2, tick=FALSE, line=yaxt.line,
      at=nY:1, labels=yaxt
    )
  }
  mtext(main, side=3, cex=par("cex.main"), font=2)
  
  # Add legend if specified
  if (plotLegend) {
    if (missing(legend.lim))
      legend.lim <- vlim
    pw <- nX + 1
    ph <- nY + 1
    addGradientLegend(
      palette, vlim, legend.lim, FALSE, legend.main,
      xlim=c(pw - 0.5 + pw*legend.line, pw - 0.5 + pw*(legend.line+0.05)), 
      ylim=c(ph/3, ph - 0.5 - ph*0.1), tick.size=legend.tick.size,
      axis.line=laxt.line, border.width=border.width
    )
  }
}

#' Plot a color palette legend
#' 
#' Add a legend to a plot window.
#' 
#' @param palette color palette.
#' @param palette.vlim limits of the values mapping to the extremities of the 
#'  color palette.
#' @param legend.vlim limits of the values to display on the legend
#' @param horizontal logical; if \code{TRUE} the legend is plotted horizontally,
#'   otherwise vertically.
#' @param main title of the legend.
#' @param xlim xlim relative to the plotting region of the rest of the plot.
#' @param ylim ylim relative to the plotting region of the rest of the plot.
#' @param tick.size size of the legend axis ticks relative to the size of the 
#'  plot window.
#' @param axis.line distance from the axis to render the axis labels as a 
#'  multiple of \code{tick.size}.
#' @param border.width line width for borders.
#' 
addGradientLegend <- function(
  palette, palette.vlim, legend.vlim, horizontal, main, xlim, ylim, 
  tick.size=0.04, axis.line=3, border.width=2
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
      legend.vlim[1], palette, palette.vlim
    )), 1),
    tail(which(palette == getColFromPalette(
      legend.vlim[2], palette, palette.vlim
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
    border="black", lwd=border.width, xpd=NA
  )
  
  # Render axis, but make sure it's balanced around 0
  labels <- c(
    seq.int(legend.vlim[1L], 0, length.out=3),
    seq.int(0, legend.vlim[2L], length.out=3)[-1]
  )
  if (length(unique(labels)) == 3)
    labels <- seq.int(legend.vlim[1L], legend.vlim[2L], length.out=5)
  labels <- prettyNum(labels, digits=2)
  if (horizontal) {
    tck <- (par("usr")[4] - par("usr")[3])*tick.size
    
    # for mapping from vlim to plot space
    v.per.x <- (xlim[2] - xlim[1])/(legend.vlim[2] - legend.vlim[1])
    zero <- xlim[1] +  v.per.x * (0 - legend.vlim[1])
    at <- c(
      seq.int(xlim[1L], zero, length.out=3),
      seq.int(zero, xlim[2L], length.out=3)[-1]
    )
    if(length(unique(at)) == 3)
      at <- seq.int(xlim[1L], xlim[2L], length.out=5)
    
    # Now plot the lines and text
    sapply(at, function(aa) {
      lines(x=c(aa, aa), y=c(ylim[1], ylim[1]-tck), lwd=border.width, xpd=NA)
    })
    text(labels, x=at, y=ylim[1]-tck*axis.line, cex=par("cex.axis"), xpd=NA)
  } else {
    tck <- (par("usr")[2] - par("usr")[1])*tick.size
    
    # for mapping from vlim to plot space
    v.per.y <- (ylim[2] - ylim[1])/(legend.vlim[2] - legend.vlim[1])
    zero <- ylim[1] +  v.per.y * (0 - legend.vlim[1])
    at <- c(
      seq.int(ylim[1L], zero, length.out=3),
      seq.int(zero, ylim[2L], length.out=3)[-1]
    )
    if(length(unique(at)) == 3)
      at <- seq.int(ylim[1L], ylim[2L], length.out=5)
    
    # draw axis ticks
    sapply(at, function(aa) {
      lines(x=c(xlim[1], xlim[1]-tck), y=c(aa, aa), lwd=border.width, xpd=NA)
    })
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
#' @param xaxt.line the number of lines into the margin at which the x axis 
#'  labels will be drawn.
#' @param maxt.line the number of lines into the margin at which the module 
#'  names will be drawn.
#' @param ylab label for the yaxis
#' @param border.width line width for borders.
#' 
plotBar <- function(
  heights, heights.lim, mas, cols, bar.width=1, drawBorders=FALSE, 
  na.col="#bdbdbd", xaxt=TRUE, plotModuleNames=TRUE, main="", xaxt.line=-0.5,
  maxt.line=3, ylab="", border.width=2
) {
  if (length(cols) == 1)
    cols <- rep(cols, length(heights))
  if (length(cols) < length(heights))
    cols <- rep(cols, ceiling(length(heights))/length(cols))[1:length(heights)]

  ylim <- heights.lim
  ylim[2] <- ylim[2] + (ylim[2] - ylim[1])*0.01
  ylim[1] <- ylim[1] - (ylim[2] - ylim[1])*0.01
  emptyPlot(
    xlim=c(0.5, length(heights)+0.5), ylim=ylim, bty="n", ylab=ylab, xpd=NA
  )
  
  # draw NAs
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
      col=cols[ii],
      border=ifelse(drawBorders, "black", NA),
      lwd=border.width
    ) 
  }
  abline(h=0, col="black", lwd=border.width)
  axis(side=2, las=2, lwd=border.width)
  
  # render module boundaries
  if (length(unique(mas)) > 1) {
    breaks <- getModuleBreaks(mas)
    for (bi in head(breaks[-1], -1)) {
      lines(x=rep(bi, 2), y=heights.lim, lwd=border.width)
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
  mtext(main, side=3, cex=par("cex.main"), font=2)
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
#' @param na.col color of missing values to plot.
#' @param yaxt logical; If \code{TRUE}, the rownames of \code{heights} will be 
#'  rendered to the left of the bars.
#' @param plotModuleNames logical; if \code{TRUE} the names of the modules are
#'  plotted along the x axis.
#' @param yaxt.line the number of lines into the margin at which the y axis 
#'  labels will be drawn.
#' @param maxt.line the number of lines into the margin at which the module 
#'  labels will be drawn.
#' @param xlab x axis label
#' @param cex.modules relative size of module names.
#' @param border.width line width for borders.
#'
plotMultiBar <- function(
  lengths, lengths.lim, cols, bar.width=1, drawBorders=FALSE, main="",
  na.col="#bdbdbd", yaxt=TRUE, plotModuleNames=TRUE, yaxt.line=0, maxt.line=2.5,
  xlab="", cex.modules=0.7, border.width=2
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
  if (length(cols) == 1) 
    cols <- matrix(cols, nrow=nrow(lengths), ncol=ncol(lengths))
  
  pw <- 0.7 # width of each plot within the 0-1 space
  
  emptyPlot(
    xlim=c(0, ncol(lengths)), ylim=c(0, nrow(lengths)*1.01), bty="n", xlab=xlab,
    xpd=NA
  )
  
  # draw NAs
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
    for (jj in seq_len(nrow(lengths))) {
      rect(
        xleft=getX(ax),
        xright=getX(lengths[jj,ii]),
        ybottom=nrow(lengths) - jj + (1 - bar.width)/2,
        ytop=nrow(lengths) - (jj - 1) - (1 - bar.width)/2,
        col=cols[jj, ii],
        border=ifelse(drawBorders, "black", NA),
        lwd=border.width
      ) 
    }
    
    # draw 0 axis
    lines(x=rep(getX(ax), 2), y=c(0, nrow(lengths)), lwd=border.width)
    
    # draw axis
    axis(
      side=1, labels=FALSE, tck=-0.025, lwd=border.width,
      at=unique(c(getX(rr[1]), getX(ax), getX(rr[2])))
    )
    axis(
      side=1, tick=FALSE, line=-0.6, las=2,
      at=unique(c(getX(rr[1]), getX(ax), getX(rr[2]))), 
      labels=prettyNum(unique(c(rr[1], ax, rr[2])), digits=2)
    )
    if (plotModuleNames) {
      mtext(
        colnames(lengths)[ii], side=3, at=ii-0.5, cex=cex.modules, font=2,
        line=maxt.line
      )
    }
  }
    
  # Draw title
  mtext(
    main, side=3, at=ncol(lengths)/2, cex=par("cex.main"), font=2, adj=0.5,
    line=1
  )
   
  # Draw sample names
  if (yaxt) {
    axis(
      side=2, tick=FALSE, las=2, at=1:nrow(lengths)-0.5,
      labels=rownames(lengths), line=yaxt.line
    )
  }
}
