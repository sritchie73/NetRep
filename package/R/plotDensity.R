#' Plot the distributions for the statistics
#'
#' @param null Three dimensional array containing the null distributions for 
#'  each of the preservation statistics for each network subset.
#' @param observed matrix containing the test statistic for each of the 
#'  preservation statistics for each network subset.
#' @importFrom scales hue_pal
#' @export 
plotStatDensities <- function(null, observed) {
  layout(cbind(matrix(1:9, 3, byrow=TRUE), 10), widths=c(0.3, 0.3, 0.3, 0.1))
  colors <- hue_pal()(nrow(null))
  for (i in seq_len(ncol(null))) {
    # Set up plot dimensions
    min_x = min(null[,i,])
    max_x = max(null[,i,])
    if (!missing(observed)) {
      min_x = min(min_x, observed[,i])
      max_x = max(max_x, observed[,i])
    }
    range = max_x - min_x
    min_x = min_x - range*0.1
    max_x = max_x + range*0.1
    ylab = ""
    xlab = ""
    par(mar=c(4, 4, 2, 1) + 0.1)
    if (i %% 3 == 1) {
      ylab="density"
    }
    if (i %in% 7:9) {
      xlab="test statistic"
    }
    for (j in seq_len(nrow(null))) {
      if (j == 1) {
        plot(0, type="n", xaxt="s", yaxt="s", main=colnames(null)[i],
             xlim=c(min_x, max_x), ylim=c(0,1.1), xlab=xlab, ylab=ylab)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
             col="#F0F0F0", border=NA)
        # Draw ggplot2 like grid
        draw.tick.lines(axTicks(1), 'x')
        draw.tick.lines(axTicks(2), 'y')
      }
      lines(density(null[j,i,]), col=colors[j])
      if (!missing(observed)) {
        abline(v=observed[j,i], col=colors[j])
      }
    }
  }
  if (i < 9) {
    while (i < 9) {
      plot.new()
      i = i + 1
    }
  }
  par(mar=c(0, 0, 2, 0) + 0.1)
  plot.new()  
  legend("top", title="subset", legend=rownames(null), fill=colors)
}

# Axis must either be 'x' or 'y'
draw.tick.lines <- function(ticks, axis) {
  dirAbline <- function(l, axis, ...) {
    switch(axis,
           "x" = abline(v=l, ...),
           "y" = abline(h=l, ...)
    )
  }
  # Major grid lines
  sapply(ticks, dirAbline, axis, col="white", lwd=2)
  # Minor grid lines
  minorlwd <- 0.5
  for(i in head(seq_along(ticks), -1)) {
    midtick <- mean(ticks[c(i, i+1)])
    dirAbline(midtick, axis, col="white", lwd=minorlwd)
  }
  # Handle boundary conditions
  lower <- switch(axis, "x" = par("usr")[1], "y" = par("usr")[2])
  upper <- switch(axis, "x" = par("usr")[3], "y" = par("usr")[4])
  difftick <- tail(ticks, 1) - midtick
  if ((ticks[1] - difftick) > lower) 
    dirAbline(ticks[1] - difftick, axis, col="white", lwd=minorlwd)
  if ((tail(ticks, 1) + midtick) < upper) 
    dirAbline(tail(ticks, 1) + difftick, axis, col="white", lwd=minorlwd)
}