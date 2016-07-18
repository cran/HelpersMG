#' symbol.Male plot a male symbol in the plotting region
#' @title Plot a male symbol in the plotting region
#' @author  Marc Girondot
#' @return Nothing
#' @param centerx The x position of the center of the circle
#' @param centery The y position of the center of the circle
#' @param rayonx The size of the rayon in the scale of the x axis
#' @param lwd The width of the line of the symbol
#' @param col The color of the symbol
#' @description Plot a male symbol in the plotting region.
#' @family Symbol
#' @examples
#' plot(x=1:2, y=c(10,20), type="n", bty="n", xlab="", ylab="")
#' 
#' rayonx <- 0.01
#' centerx <- 1.2
#' centery <- 15
#' 
#' symbol.Male(centerx=centerx, centery = centery, rayonx=rayonx)
#' symbol.Female(centerx=centerx+0.5, centery = centery, rayonx=rayonx)
#' 
#' rayonx <- 0.03
#' centerx <- 1.2
#' centery <- 18
#' 
#' symbol.Male(centerx=centerx, centery = centery, rayonx=rayonx, lwd=3)
#' symbol.Female(centerx=centerx+0.5, centery = centery, rayonx=rayonx, lwd=3, col="red")
#' 
#' rayonx <- 0.05
#' centerx <- 1.4
#' centery <- 13
#' 
#' symbol.Male(centerx=centerx, centery = centery, rayonx=rayonx, lwd=4, col="blue")
#' symbol.Female(centerx=centerx+0.5, centery = centery, rayonx=rayonx, lwd=4, col="red")
#' @export


symbol.Male <- function(centerx, centery, rayonx, lwd=2, col="black") {
xr <- ScalePreviousPlot()$xlim["range"]
yr <- ScalePreviousPlot()$ylim["range"]
ratio <- par("pin")[1]/par("pin")[2]
rayony <- rayonx*(yr/xr)*ratio

angle <- seq(from=0, to=2*pi, length.out = 20)
x <- centerx+rayonx*cos(angle)
y <- centery+rayony*sin(angle)
segments(x0=x, y0=y, x1=x[c(20, 1:19)], y1=y[c(20, 1:19)], col=col, lwd=lwd)

x0 = centerx+rayonx*sqrt(2)/2
y0 = centery+rayony*sqrt(2)/2
x1 = centerx+3*rayonx*sqrt(2)/2
y1 = centery+3*rayony*sqrt(2)/2
segments(x0=x0, y0=y0, x1=x1, y1=y1, col=col, lwd=lwd)

x11 <- centerx+2*rayonx*cos(2*pi*60/360)
y11 <- centery+2*rayony*sin(2*pi*60/360)

x12 <- centerx+2.5*rayonx*cos(2*pi*45/360)
y12 <- centery+2.5*rayony*sin(2*pi*45/360)

x13 <- centerx+2*rayonx*cos(2*pi*30/360)
y13 <- centery+2*rayony*sin(2*pi*30/360)


x <- c(x1, x11, x12, x13, x1)
y <- c(y1, y11, y12, y13, y1)
polygon(x, y, col=col, lwd=lwd, border=col)
}