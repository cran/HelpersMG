#' plot_errbar plot a xy graph with error bar on x and/or y
#' @title Plot a xy graph with error bar on x and/or y
#' @author Marc Girondot
#' @return Nothing
#' @param ... Parameters for plot() such as main= or ylim=
#' @param errbar.x The length of error bars for x. Recycled if necessary.
#' @param errbar.x.plus The length of positive error bars for x. Recycled if necessary.
#' @param errbar.x.minus The length of negative error bars for x. Recycled if necessary.
#' @param errbar.y The length of error bars for y. Recycled if necessary.
#' @param errbar.y.plus The length of positive error bars for y. Recycled if necessary.
#' @param errbar.y.minus The length of negative error bars for y. Recycled if necessary.
#' @param x.plus The absolut position of the positive error bar for x. Recycled if necessary.
#' @param x.minus The absolut position of the negative error bar for x. Recycled if necessary.
#' @param y.plus The absolut position of the positive error bar for y. Recycled if necessary.
#' @param y.minus The absolut position of the nagative error bar for y. Recycled if necessary.
#' @param errbar.tick Size of small ticks at the end of error bars defined as a proportion of total width or height graph size.
#' @param errbar.lwd Error bar line width, see par("lwd")
#' @param errbar.lty Error bar line type, see par("lwd")
#' @param errbar.col Error bar line color, see par("col")
#' @param errbar.y.polygon If true, the errors are shown as a filed polygon.
#' @param errbar.y.polygon.list List of parameters to be used for polygon.
#' @param names The names of the points to be used with show_name().
#' @param add If true, add the graph to the previous one.
#' @family plot and barplot functions
#' @seealso \code{barplot_errorbar}
#' @description To plot data, just use it as a normal plot but add the errbar.x 
#' and errbar.y values or errbar.x.minus, errbar.x.plus if bars for x axis are 
#' asymetric and errbar.y.minus, errbar.y.plus if bars for y axis are 
#' asymetric. Use x.plus, x.minus, y.plus and y.minus to set absolut limits for
#' error bars. Note that x.plus and x.minus have priority over errbar.x, errbar.x.minus and
#' errbar.x.plus and that y.plus and y.minus have priority over errbar.y, errbar.y.minus and
#' errbar.y.plus.\cr
#' The parameter errbar.y.polygon=TRUE permits to define error as an envolop for y axis.
#' @examples
#' \dontrun{
#' plot_errbar(1:100, rnorm(100, 1, 2), 
#'		xlab="axe x", ylab="axe y", bty="n", xlim=c(1,100), 
#' 		errbar.x=2, errbar.y=rnorm(100, 1, 0.1))
#' x <- 1:100
#' plot_errbar(x=1:100, rnorm(100, 1, 2), 
#'                	xlab="axe x", ylab="axe y", bty="n", xlim=c(1,100), 
#'             		x.minus=x-2, x.plus=x+2)
#' x <- 1:100
#' plot_errbar(x=1:100, rnorm(100, 1, 2), 
#'                	xlab="axe x", ylab="axe y", bty="n", 
#'                	pch=21, bg="white", 
#'             		x.minus=x-10, x.plus=x+10)
#' x <- (1:200)/10
#' y <- sin(x)
#' plot_errbar(x=x, y=y, xlab="axe x", ylab="axe y", bty="n", xlim=c(1,20), 
#'      y.minus=y-1, y.plus=y+1, ylim=c(-3, 3), type="l",  
#' 		errbar.y.polygon=TRUE, 
#' 		errbar.y.polygon.list=list(border=NA, col=rgb(0, 0, 0, 0.5)))
#' 		}
#' @export


plot_errbar <- function(..., 
                        errbar.x=NULL, errbar.y=NULL, 
                        errbar.x.plus=NULL, errbar.x.minus=NULL, 
                        errbar.y.plus=NULL, errbar.y.minus=NULL,
                        x.plus=NULL, x.minus=NULL,
                        y.plus=NULL, y.minus=NULL,
                        errbar.tick=1/50, 
                        errbar.lwd=par("lwd"), 
                        errbar.lty=par("lty"), 
                        errbar.col=par("fg"), 
                        errbar.y.polygon=FALSE, 
                        errbar.y.polygon.list=list(NULL), 
                        names=NULL, 
                        add=FALSE) 
{
  
  # errbar.x=NULL; errbar.y=NULL; errbar.x.plus=NULL; errbar.x.minus=NULL; errbar.y.plus=NULL; errbar.y.minus=NULL; x.plus=NULL; x.minus=NULL; y.plus=NULL; y.minus=NULL; errbar.tick=1/50; errbar.lwd=par("lwd"); errbar.lty=par("lty"); errbar.col=par("fg"); errbar.y.polygon=FALSE; errbar.y.polygon.list=list(NULL); add=FALSE
  # par.plot <- list(x=x.axis, y=CTE, las=1, type="l", xlim=c(as.Date("1997-01-01"), as.Date("2014-01-01")), ylim=c(28, 32), bty="n", xlab="Year", ylab=expression("constant incubation temperature (" *degree*"C)"), xaxt="n")
  # y.plus=CTE.plus
  # y.minus=CTE.moins
  # errbar.y.polygon=TRUE
  # errbar.y.polygon.list=list(border=NA, col=rgb(0, 0, 0, 0.5))
  # par.plot <- list(1:100, rnorm(100, 1, 2),xlab="axe x", ylab="axe y", bty="n", xlim=c(1,100));errbar.x=2; errbar.y=rnorm(100, 1, 0.1)
  
  par.plot <- list(...)
  
  par(xpd=FALSE)
  
  x <- par.plot[["x"]]
  if (is.null(x)) {
    x <- par.plot[[1]]
    names(par.plot)[1] <- "x"
  }
  
  if (is.data.frame(x) | is.matrix(x)) {
    y <- x[,2]
    x <- x[,1]
  } else {
    y <- par.plot[["y"]]
  }
  if (is.null(y)) {
    y <- par.plot[[2]]
    names(par.plot)[2] <- "y"
  }
  
  if (!is.null(x.plus)) errbar.x.plus <- x.plus-x
  if (!is.null(x.minus)) errbar.x.minus <- x-x.minus
  if (!is.null(y.plus)) errbar.y.plus <- y.plus-y
  if (!is.null(y.minus)) errbar.y.minus <- y-y.minus
  
  if (is.null(errbar.x.minus) & !is.null(errbar.x)) {
    errbar.x.minus <- errbar.x
  }
  if (is.null(errbar.x.plus) & !is.null(errbar.x)) {
    errbar.x.plus <- errbar.x
  }
  if (is.null(errbar.y.minus) & !is.null(errbar.y)) {
    errbar.y.minus <- errbar.y
  }
  if (is.null(errbar.y.plus) & !is.null(errbar.y)) {
    errbar.y.plus <- errbar.y
  }
  
  if (add) {
    # Si je superpose le graphique à un précédent:
    s <- ScalePreviousPlot()
    par(new=TRUE)
    # Je retire les axes et je fixe xlim et ylim
    pp <- modifyList(par.plot, list(xlim=s$xlim[1:2], ylim=s$ylim[1:2], xlab="", ylab="", main="", axes=FALSE))
    do.call(plot, modifyList(pp, list(type="n")))
    
  } else {
    # Je fais un nouveau graphique
    # Là il n'y a pas de xlim, ylim; je le laisse calculer
    pp <- par.plot
    pp <- modifyList(pp, list(x=c(min(x-ifelse(is.null(errbar.x.minus), 0, errbar.x.minus)), 
                                  max(x+ifelse(is.null(errbar.x.plus), 0, errbar.x.plus)))
    ))
    pp <- modifyList(pp, list(y=c(min(y-ifelse(is.null(errbar.y.minus), 0, errbar.y.minus)), 
                                  max(y+ifelse(is.null(errbar.y.plus), 0, errbar.y.plus)))
    ))
    
    do.call(plot, modifyList(pp, list(type="n")))
    
    # Je rajoute xlim et ylim sur par.plot: 2020-04-02
    
    s <- ScalePreviousPlot()
    par.plot <- modifyList(par.plot, list(xlim=s$xlim[1:2], ylim=s$ylim[1:2]))
  }
  
  
  
  
  if (errbar.y.polygon) {
    # je dois faire un polygon
    # Dans ce cas, c'est bon
    vx <- c(x, rev(x))
    vy <- c(y-errbar.y.minus, rev(y+errbar.y.plus))
    errbar.y.polygon.list <- modifyList(errbar.y.polygon.list, list(x=vx, y=vy))
    do.call(polygon, errbar.y.polygon.list)
    
    
  } else {
    
    sizebar <- (par("usr")[4]-par("usr")[3])*errbar.tick
    
    # Je fais les barres d'erreur
    if (!is.null(errbar.x.minus)) {
      segments(x-errbar.x.minus, y, x, y, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
      segments(x-errbar.x.minus, y-sizebar, x-errbar.x.minus, y+sizebar, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    }
    if (!is.null(errbar.x.plus)) {
      segments(x+errbar.x.plus, y, x, y, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
      segments(x+errbar.x.plus, y-sizebar, x+errbar.x.plus, y+sizebar, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    }
    
    sizebar <- (par("usr")[2]-par("usr")[1])*errbar.tick
    
    if (!is.null(errbar.y.minus)) {
      segments(x, y-errbar.y.minus, x, y, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
      segments(x-sizebar, y-errbar.y.minus, x+sizebar, y-errbar.y.minus, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    }
    if (!is.null(errbar.y.plus)) {
      segments(x, y+errbar.y.plus, x, y, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
      segments(x-sizebar, y+errbar.y.plus, x+sizebar, y+errbar.y.plus, 
               col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    }
  }
  
  # Et là je fais le plot
  do.call(plot_add, par.plot)
  
  if (is.null(names)) names <- paste0(as.character(x), ";", as.character(y))
  
  return(invisible(list(x=x, y=y, names=names)))
  
}
