#' Plot a symmetric heatmap as a triangle
#' 
#' @param values values to plot on the heatmap
#' @param palette color palette to interpolate over
#' @param vlim range of values to use when mapping values to the \code{palette}.
#' @param mas ordered subset of the moduleAssignments vector
#' @param na.indices indices of missing values to plot
#' @param na.col color of missing values to plot.
#' 
plotTriangleHeatmap <- function(
  values, palette, vlim, mas, na.indices=NULL, na.col="#bdbdbd"
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
        y=c(0, 0, height, 0)
      )
    }
  }
  
  # render border of plot
  polygon(
    x=c(0.5, nGenes+0.5, nGenes/2+0.5, 0.5), y=c(0, 0, nGenes/2, 0),
    lwd=2, xpd=TRUE
  )
}

#' Plot a square heatmap
#' 
#' @param values values to plot on the heatmap
#' @param palette color palette to interpolate over
#' @param vlim range of values to use when mapping values to the \code{palette}.
#' @param mas ordered subset of the moduleAssignments vector
#' @param na.indices indices of missing values to plot
#' @param na.col color of missing values to plot.
#' 
plotSquareHeatmap <- function(
  values, palette, vlim, mas, na.indices=NULL, na.col="#bdbdbd"
) {
  nGenes <- ncol(values) + length(na.indices)
  emptyPlot(xlim=c(0.5, nGenes+0.5), ylim=c(0.5, nGenes+0.5), bty="n")
  palette <- colorRampPalette(palette)(255)
  
  # render squares / triangles
  ci <- 1
  for (ii in 1:nGenes) {
    cj <- 1
    for (jj in 1:nGenes) {
      if (ii %nin% na.indices && jj %nin% na.indices) {
        col <- getColFromPalette(values[ci, cj], palette, vlim)
        cj <- cj + 1
      } else {
        col <- na.col
      }
      rect(
        xleft = jj - 0.5,
        xright = jj + 0.5,
        ybottom = (nGenes - (ii - 1)) - 0.5,
        ytop = (nGenes - (ii - 1)) + 0.5,
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
    for (mi in seq_along(breaks)[-1]) {
      rect(
        xleft = breaks[mi - 1],
        xright = breaks[mi],
        ybottom = (nGenes + 0.5) - (breaks[mi] - 0.5),
        ytop = (nGenes + 0.5) - (breaks[mi - 1] - 0.5),
        border="black"
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
    xpd=TRUE,
    lwd=2
  )
}

#' Plot a color palette legend
#' 
#' Legend will fill an entire plot window.
#' 
#' @param palette color palette.
#' @param palette.vlim limits of the values mapping to the extremities of the 
#'  color palette.
#' @param legend.vlim limits of the values to display on the legend
#' @param horizontal logical; if \code{TRUE} the legend is plotted horizontally,
#'   otherwise vertically.
#' @param main title of the legend.
#' 
plotGradientLegend <- function(
  palette, palette.vlim, legend.vlim, horizontal, main
) {

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
#' @param drawBorder logical; if \code{TRUE} a border is drawn around each bar.
#' @param horizontal logical; if \code{TRUE} bars are plotted horizontally.
#' 
plotBar <- function(
  heights, heights.lim, mas, cols, bar.width=1, drawBorder=FALSE, 
  horizontal=FALSE
) {
  if (length(cols) == 1) 
    cols <- rep(cols, length(heights))
  if (horizontal) {
    emptyPlot(xlim=heights.lim, ylim=c(0.5, length(heights)+0.5), bty="n")
    for (ii in 1:length(heights)) {
      rect(
        xleft=0,
        xright=heights[ii],
        ybottom=length(heights)-(ii-1)-bar.width/2,
        ytop=length(heights)-(ii-1)+bar.width/2,
        col=cols[ii],
        border=ifelse(drawBorder, "black", NA),
        lwd=2
      )
    }
    abline(v=0, col="black", lwd=2)
    axis(side=1, lwd=2)
  } else {
    emptyPlot(xlim=c(0.5, length(heights)+0.5), ylim=heights.lim, bty="n")
    for (ii in 1:length(heights)) {
      rect(
        xleft=ii-bar.width/2,
        xright=ii+bar.width/2,
        ybottom=0,
        ytop=heights[ii],
        col=cols[ii],
        border=ifelse(drawBorder, "black", NA),
        lwd=2
      ) 
    }
    abline(h=0, col="black", lwd=2)
    axis(side=2, las=2, lwd=2)
  }
  
  # render module boundaries
  if (length(unique(mas)) > 1) {
    breaks <- getModuleBreaks(mas)
    if (horizontal) {
      abline(h=length(heights) + 0.5 - breaks, lwd=2)
    } else {
      abline(v=breaks, lwd=2)
    }
  }
}
