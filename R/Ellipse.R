#' ellipse plots an ellipse
#' @title Plot an ellipse
#' @author marc.girondot@@u-psud.fr
#' @return Nothing
#' @param center.x Center of the ellipse on x axis
#' @param center.y Center of the ellipse on y axis
#' @param radius.x Radius along the x axis
#' @param radius.y Radius along the y axis
#' @param radius.x.lower Radius along the x axis, at left of center
#' @param radius.x.upper Radius along the x axis, at right of center
#' @param radius.y.lower Radius along the y axis, at bottom of center
#' @param radius.y.upper Radius along the y axis, at top of center
#' @param alpha Rotation in radians
#' @param binconf.x A data.frame or a matrix with two columns, x and n or with three columns, PointEst, Lower, and Upper
#' @param binconf.y A data.frame or a matrix with two columns, x and n or with three columns, PointEst, Lower, and Upper
#' @param control.binconf A list with options for binomial confidence
#' @param length Number of points to draw the ellipse
#' @param ... Graphical parameters
#' @description Plot a ellipse dined by the center and the radius. The options for 
#' binomial confidence are:\cr
#' - alpha is 1 - confidence interval\cr
#' - method  must be one of these "wilson", "exact", "asymptotic"\cr
#' col parameter can be a list of colors. See examples
#' @examples
#' plot(0:1, 0:1, xlim=c(0, 1), ylim=c(0,1), lty=2, type="l", las=1, bty="n", 
#'      xlab="Variable x", ylab="variable y")
#'  
#' ellipse(center.x = c(0.2, 0.3, 0.25), center.y = c(0.7, 0.6, 0.55), 
#'         radius.x = c(0.1, 0.1, 0.1), radius.y = c(0.15, 0.2, 0.4), 
#'         border=NA, col=rgb(red = 0.1, green = 0.1, blue = 0.1, alpha = 0.1))
#' 
#' ellipse(center.x = 0.5, center.y = 0.5, 
#'         radius.x.lower = 0.1, radius.x.upper = 0.3, 
#'         radius.y = 0.2, 
#'         border=NA, col=rgb(red = 0.1, green = 0.1, blue = 0.1, alpha = 0.1))
#' 
#' ellipse(center.x = 0.6, center.y = 0.3, 
#'         radius.x.lower = 0.3, radius.x.upper = 0.3, 
#'         radius.y.lower = 0.2, radius.y.upper = 0.4, 
#'         border=NA, col=rgb(red = 0.1, green = 0.1, blue = 0.1, alpha = 0.1))
#' 
#' plot(0:1, 0:1, xlim=c(0, 1), ylim=c(0,1), lty=2, type="l", bty="n", asp=1, 
#'      xlab="Variable x", ylab="variable y", axes=FALSE)
#' axis(1, at=c(0, 0.25, 0.5, 0.75, 1))
#' axis(2, at=c(0, 0.25, 0.5, 0.75, 1), las=1)
#' 
#' ellipse(center.x = 0.5, center.y = 0.5, radius.x = 0.2, radius.y = 0.4, 
#'        border=NA, col=rgb(red = 0.1, green = 0.1, blue = 0.1, alpha = 0.1))
#' ellipse(center.x = 0.5, center.y = 0.5, radius.x = 0.2, radius.y = 0.4, 
#'         border=NA, col=rgb(red = 0.1, green = 0.1, blue = 0.1, alpha = 0.1), alpha = pi/4)
#' 
#' plot(0:1, 0:1, xlim=c(0, 1), ylim=c(0,1), lty=2, type="l", las=1, bty="n", 
#'      xlab="Variable x", ylab="variable y")
#' 
#' for (k in 0:8)
#'   ellipse(center.x=0.5, center.y=0.5, radius.x=0.1, radius.y=0.4, 
#'           alpha=seq(from=0, to=pi/4, length=9)[k], 
#'           border=rainbow(9)[k])
#'
#' # Exemple with confidence of proportions
#' males <- c(10, 25, 3, 4)
#' N <- c(12, 52, 17, 10)
#' 
#' males2 <- c(12, 20, 3, 6)
#' N2 <- c(15, 50, 20, 12)
#' 
#' plot(0:1, 0:1, xlim=c(0, 1), ylim=c(0,1), lty=2, type="l", las=1, bty="n", 
#'      xlab="Variable x", ylab="variable y")
#'
#' ellipse(binconf.x = data.frame(x=males, n=N), binconf.y = data.frame(x=males2, n=N2),  
#'         border=NA, col=rgb(red = 0.1, green = 0.5, blue = 0.1, alpha = 0.1))
#'         
#' plot(0:1, 0:1, xlim=c(0, 1), ylim=c(0,1), lty=2, type="l", las=1, bty="n", 
#'      xlab="Variable x", ylab="variable y")
#'      
#' ellipse(binconf.x = data.frame(x=males, n=N), 
#'         binconf.y = data.frame(PointEst=c(0.1, 0.2, 0.3, 0.5), 
#'                                Lower=c(0.02, 0.12, 0.25, 0.30), 
#'                                Upper=c(0.18, 0.29, 0.35, 0.67)), 
#'         border=NA, col=rgb(red = 0.1, green = 0.5, blue = 0.1, alpha = 0.1))
#'         
#' # Examples with a gradient
#' plot(0:1, 0:1, xlim=c(0, 1), ylim=c(0,1), lty=2, type="l", las=1, bty="n", 
#'      xlab="Variable x", ylab="variable y")
#' ellipse(center.x = 0.6, center.y = 0.3, 
#'         radius.x.lower = 0.3, radius.x.upper = 0.3, 
#'         radius.y.lower = 0.2, radius.y.upper = 0.4, 
#'         border=NA, col=grey.colors(100, alpha = 0.1))
#'         
#' plot(0:1, 0:1, xlim=c(0, 1), ylim=c(0,1), lty=2, type="l", las=1, bty="n", 
#'      xlab="Variable x", ylab="variable y")
#' ellipse(binconf.x = data.frame(x=males, n=N), binconf.y = data.frame(x=males2, n=N2),  
#'         border=NA, col=grey.colors(100, alpha = 0.1))
#' 
#' @export


ellipse <- function(center.x = 0, center.y = 0, 
                        radius.x = 1, radius.y = 1, 
                        radius.x.lower=NULL, radius.x.upper=NULL, 
                        radius.y.lower=NULL, radius.y.upper=NULL, 
                    alpha = 0, 
                    binconf.x=NULL, 
                    binconf.y=NULL, 
                    control.binconf=list(alpha = 0.05, method = "wilson"), 
                        length=100, ...) {
  
  p3p <- list(...)
  
  if (!is.null(binconf.x)) {
    if ((class(binconf.x)=="binconf") | (ncol(binconf.x) == 3)) {
      bc.x <- binconf.x
    } else {
      control.binconf.x <- modifyList(c(control.binconf, list(x=binconf.x[, "x"]), 
                                        list(n=binconf.x[, "n"])), 
                                      list(include.x = FALSE, include.n = FALSE, return.df = FALSE))
      bc.x <- do.call(getFromNamespace(".BinomialConfidence", ns="HelpersMG"), control.binconf.x) 
    }
    
    center.x = bc.x[, "PointEst"]
    radius.x.lower = bc.x[, "PointEst"]-bc.x[, "Lower"]
    radius.x.upper = bc.x[, "Upper"]-bc.x[, "PointEst"] 
  }
  
  if (!is.null(binconf.y)) {
    if ((class(binconf.y)=="binconf") | (ncol(binconf.y) == 3))  {
      bc.y <- binconf.y
    } else {
      control.binconf.y <- modifyList(c(control.binconf, list(x=binconf.y[, "x"]), 
                                        list(n=binconf.y[, "n"])), 
                                      list(include.x = FALSE, include.n = FALSE, return.df = FALSE))
      bc.y <- do.call(getFromNamespace(".BinomialConfidence", ns="HelpersMG"), control.binconf.y) 
    }
    
    center.y = bc.y[, "PointEst"]
    radius.y.lower = bc.y[, "PointEst"]-bc.y[, "Lower"] 
    radius.y.upper = bc.y[, "Upper"]-bc.y[, "PointEst"] 
  }
  
  if (is.null(radius.x.lower))  radius.x.lower <- radius.x
  if (is.null(radius.y.lower))  radius.y.lower <- radius.y
  if (is.null(radius.x.upper))  radius.x.upper <- radius.x
  if (is.null(radius.y.upper))  radius.y.upper <- radius.y
  
  ncol <- 1
  if (!is.null(p3p$col)) {
      ncol <- length(p3p$col)
      p3p$col <- rev(p3p$col)
  }
  
  
  for (k in seq_along(center.x)) {
    
    cptcol <- 1
    for (nc in seq(from=1, to=1/ncol, length.out = ncol)) {
 
  theta <- seq(0, pi / 2, length=length/4)
  x <- center.x[k] + radius.x.upper[k] * nc * cos(theta) * cos(alpha) - radius.y.upper[k] * nc * sin(theta) * sin(alpha)
  y <- center.y[k] + radius.x.upper[k] * nc * cos(theta) * sin(alpha) + radius.y.upper[k] * nc * sin(theta) * cos(alpha)
  
  theta <- seq(pi / 2, pi, length=length/4)
  x <- c(x, center.x[k] + radius.x.lower[k] * nc * cos(theta) * cos(alpha) - radius.y.upper[k] * nc * sin(theta) * sin(alpha))
  y <- c(y, center.y[k] + radius.x.lower[k] * nc * cos(theta) * sin(alpha) + radius.y.upper[k] * nc * sin(theta) * cos(alpha))
  
  theta <- seq(pi, 3/2*pi, length=length/4)
  x <- c(x, center.x[k] + radius.x.lower[k] * nc * cos(theta) * cos(alpha) - radius.y.lower[k] * nc * sin(theta) * sin(alpha))
  y <- c(y, center.y[k] + radius.x.lower[k] * nc * cos(theta) * sin(alpha) + radius.y.lower[k] * nc * sin(theta) * cos(alpha))
  
  theta <- seq(3/2*pi, 2 * pi, length=length/4)
  x <- c(x, center.x[k] + radius.x.upper[k] * nc * cos(theta) * cos(alpha) - radius.y.lower[k] * nc * sin(theta) * sin(alpha))
  y <- c(y, center.y[k] + radius.x.upper[k] * nc * cos(theta) * sin(alpha) + radius.y.lower[k] * nc * sin(theta) * cos(alpha))
  
  if ((ncol == 1)) {
    do.call(polygon, modifyList(p3p, list(x=x, y=y)))
  } else {
    p3p_ec <- modifyList(p3p, list(col=p3p$col[cptcol]))
    do.call(polygon, modifyList(p3p_ec, list(x=x, y=y)))
  }
  cptcol <- cptcol + 1
    }
 
  }
}

