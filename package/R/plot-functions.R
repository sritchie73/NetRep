#' Plot a symmetric heatmap as a triangle
#' 
#' @param values values to plot on the heatmap
#' @param palette color palette to interpolate over
#' @param vlim range of values to use when mapping values to the \code{palette}.
#' @param mas ordered subset of the moduleAssignments vector
#' 
plotTriangleHeatmap <- function(values, palette, vlim, mas) {
  nGenes <- ncol(values)
  emptyPlot(xlim=c(0.5, nGenes + 0.5), ylim=c(0, nGenes/2), bty="n")
  palette <- colorRampPalette(palette)(255)
  
  # render squares / triangles
  for (ii in 1:nGenes) {
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
      
      col <- getColFromPalette(values[ii, jj], palette, vlim)
      polygon(
        x=c(leftx, leftx+0.5, rightx, leftx+0.5, leftx),
        y=c(topy-0.5, topy, topy-0.5, boty, topy-0.5),
        col=col, border=col
      )
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
    for (mm in 2:length(breaks)) {
      mheight <- (breaks[mm] - breaks[mm - 1])/2
      halfway <-  mheight + breaks[mm - 1]
    }
  }
  
  # render border of plot
  polygon(
    x=c(0.5, nGenes+0.5, nGenes/2+0.5, 0.5), y=c(0, 0, nGenes/2, 0),
    lwd=2, xpd=TRUE
  )
}
