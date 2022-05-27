#' Show the name of a point
#' @title Show the name of a point
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Name of the point
#' @param points A list with x, y and names elements
#' @param x The x coordinates
#' @param y The y coordinates.
#' @param names The names of the points
#' @param col Color of the legend.
#' @param silent TRUE or FALSE
#' @family plot and barplot functions
#' @seealso \code{plot_errorbar}
#' @description Click on a point in plot region and it will tell you what is the point.
#' @examples
#' \dontrun{
#' k <- plot_errbar(1:100, rnorm(100, 1, 2), 
#'		xlab="axe x", ylab="axe y", bty="n", xlim=c(1,100), 
#' 		errbar.x=2, errbar.y=rnorm(100, 1, 0.1))
#' show_name(k)
#' k <- plot_errbar(1:10, rnorm(10, 1, 2), 
#'		xlab="axe x", ylab="axe y", bty="n", xlim=c(1,10), 
#' 		errbar.x=2, errbar.y=rnorm(10, 1, 0.1), 
#' 		names=LETTERS[1:10])
#' show_name(k)
#' k <- plot_errbar(1:10, rnorm(10, 1, 2), 
#'		xlab="axe x", ylab="axe y", bty="n", xlim=c(1,10), 
#' 		errbar.x=2, errbar.y=rnorm(10, 1, 0.1))
#' show_name(k, names=LETTERS[1:10])
#' 		}
#' @export


show_name <- function(points=NULL, x=NULL, y=NULL, names=NULL, col="red", 
                      silent=FALSE) {
  
  if (!silent) message("Click on a point in the plot")
  if (!is.null(points)) {
    if (is.null(x)) x <- points[["x"]]
    if (is.null(y)) y <- points[["y"]]
    if (is.null(names)) names <- points[["names"]]
  }
  
  k <- locator(n=1)
  kpos <- which.min(abs(x - k$x) + 
                      abs(y - k$y))
  
  dx <- ScalePreviousPlot()$x["range"]/30
  dy <- ScalePreviousPlot()$y["range"]/30
  
  if (x[kpos] < ScalePreviousPlot()$x["center"]) {
    if (y[kpos] < ScalePreviousPlot()$y["center"]) {
      par(xpd=TRUE)
      text(x=x[kpos]+dx, 
           y=y[kpos]+dy, 
           labels = names[kpos], 
           col=col, 
           pos=4)
      segments(x0=x[kpos], x1=x[kpos]+dx, 
               y0=y[kpos], y1=y[kpos]+dy, col=col)
    } else {
      par(xpd=TRUE)
      text(x=x[kpos]+dx, 
           y=y[kpos]-dy, 
           labels = names[kpos], 
           col=col, 
           pos=4)
      segments(x0=x[kpos], x1=x[kpos]+dx, 
               y0=y[kpos], y1=y[kpos]-dy, col=col)
    }
  } else {
    if (y[kpos] < ScalePreviousPlot()$y["center"]) {
      par(xpd=TRUE)
      text(x=x[kpos]-dx, 
           y=y[kpos]+dy, 
           labels = names[kpos], 
           col=col, 
           pos=2)
      segments(x0=x[kpos], x1=x[kpos]-dx, 
               y0=y[kpos], y1=y[kpos]+dy, col=col)
    } else {
      par(xpd=TRUE)
      text(x=x[kpos]-dx, 
           y=y[kpos]-dy, 
           labels = names[kpos], 
           col=col, 
           pos=2)
      segments(x0=x[kpos], x1=x[kpos]-dx, 
               y0=y[kpos], y1=y[kpos]-dy, col=col)
    }
  }
  return(invisible(names[kpos]))
  
}
