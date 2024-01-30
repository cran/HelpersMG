#' cArrows draws curved lines with arrowhead
#' @title Draw curved lines with arrowhead
#' @author Modified from iGraph
#' @return A list wit lab.x and lab.y being the position where to draw label
#' @param x1 coordinates of points from which to draw.
#' @param y1 coordinates of points from which to draw.
#' @param x2 coordinates of points to which to draw.
#' @param y2 coordinates of points to which to draw.
#' @param code integer code (1, 2, or 3), determining kind of arrows to be drawn.
#' @param size size of the arrowhead.
#' @param width width of the arrowhead.
#' @param open shape of the arrowhead.
#' @param sh.adj Shift the beginning of the line.
#' @param sh.lwd width of the line.
#' @param sh.col color of the line.
#' @param sh.lty type of line.
#' @param h.col color of the arrowhead.
#' @param h.col.bo color of the arrowhead border.
#' @param h.lwd width of the arrowhead.
#' @param h.lty type of line for the arrowhead.
#' @param curved 0 is a straigth line, positive of negative value make the line curved.
#' @param beautiful.arrow if open is false, make the arrowhead more beautiful.
#' @description Draw a curved line with arrowhead.
#' @examples
#' plot(c(1, 10), c(1, 10), type="n", bty="n")
#' cArrows(x1=2, y1=2, x2=6, y2=6, curved=1)
#' cArrows(x1=2, y1=2, x2=6, y2=6, curved=0)
#' cArrows(x1=2, y1=2, x2=6, y2=6, curved=1, sh.adj=1)
#' cArrows(x1=2, y1=2, x2=6, y2=6, curved=-1, open=FALSE)
#' cArrows(x1=9, y1=2, x2=6, y2=6, curved=-1, open=FALSE, sh.col="red")
#' cArrows(x1=9, y1=9, x2=6, y2=6, curved=-1, open=FALSE, h.col="red")
#' cArrows(x1=2, y1=9, x2=6, y2=6, curved=1, open=FALSE, h.col="red", h.col.bo="red")
#' @export


cArrows <- function (x1, y1, x2, y2, code = 2, size = 1, width = 1.2/4/cin, 
                     open = TRUE, sh.adj = 0.1, sh.lwd = 1, sh.col = par("fg"), 
                     sh.lty = 1, h.col = sh.col, h.col.bo = sh.col, h.lwd = sh.lwd, 
                     h.lty = sh.lty, curved = FALSE, beautiful.arrow=2/3) 
{
  
  # plot(c(1, 10), c(1, 10))  
  #  x1=1; y1=1; x2=5; y2=5; code=2; size=1; width=1.2/4/0.2; open=TRUE; sh.adj = 0.1; sh.lwd = 1; sh.col = par("fg"); sh.lty = 1; h.col = sh.col; h.col.bo = sh.col; h.lwd = sh.lwd; h.lty = sh.lty; curved = FALSE
  
  
  cin <- size * par("cin")[2]
  width <- width * (1.2/4/cin)
  uin <- 1/xyinch()
  x <- sqrt(seq(0, cin^2, length = floor(35 * cin) + 2))
  delta <- sqrt(h.lwd) * par("cin")[2] * 0.005
  x.arr <- c(-rev(x), -x)
  wx2 <- width * x^2
  y.arr <- c(-rev(wx2 + delta), wx2 + delta)
  # Pourquoi rajouter NA ?
  deg.arr <- c(atan2(y.arr, x.arr), NA)
  r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)
  bx1 <- x1
  bx2 <- x2
  by1 <- y1
  by2 <- y2
  lx <- length(x1)
  r.seg <- rep(cin * sh.adj, lx)
  theta1 <- atan2((y1 - y2) * uin[2], (x1 - x2) * uin[1])
  th.seg1 <- theta1 + rep(atan2(0, -cin), lx)
  theta2 <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
  th.seg2 <- theta2 + rep(atan2(0, -cin), lx)
  x1d <- y1d <- x2d <- y2d <- 0
  if (code %in% c(1, 3)) {
    x2d <- r.seg * cos(th.seg2)/uin[1]
    y2d <- r.seg * sin(th.seg2)/uin[2]
  }
  if (code %in% c(2, 3)) {
    x1d <- r.seg * cos(th.seg1)/uin[1]
    y1d <- r.seg * sin(th.seg1)/uin[2]
  }
  if (is.logical(curved) && all(!curved) || is.numeric(curved) && 
      all(!curved)) {
    segments(x1 + x1d, y1 + y1d, x2 + x2d, y2 + y2d, lwd = sh.lwd, 
             col = sh.col, lty = sh.lty)
    phi <- atan2(y1 - y2, x1 - x2)
    r <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
    lc.x <- x2 + 2/3 * r * cos(phi)
    lc.y <- y2 + 2/3 * r * sin(phi)
  }  else {
    if (is.numeric(curved)) {
      lambda <- curved
    } else {
      lambda <- as.logical(curved) * 0.5
    }
    lambda <- rep(lambda, length.out = length(x1))
    c.x1 <- x1 + x1d
    c.y1 <- y1 + y1d
    c.x2 <- x2 + x2d
    c.y2 <- y2 + y2d
    midx <- (x1 + x2)/2
    midy <- (y1 + y2)/2
    spx <- midx - lambda * 1/2 * (c.y2 - c.y1)
    spy <- midy + lambda * 1/2 * (c.x2 - c.x1)
    sh.col <- rep(sh.col, length = length(c.x1))
    sh.lty <- rep(sh.lty, length = length(c.x1))
    sh.lwd <- rep(sh.lwd, length = length(c.x1))
    lc.x <- lc.y <- numeric(length(c.x1))
    for (i in seq_len(length(c.x1))) {
      if (lambda[i] == 0) {
        segments(c.x1[i], c.y1[i], c.x2[i], c.y2[i], 
                 lwd = sh.lwd[i], col = sh.col[i], lty = sh.lty[i])
        phi <- atan2(y1[i] - y2[i], x1[i] - x2[i])
        r <- sqrt((x1[i] - x2[i])^2 + (y1[i] - y2[i])^2)
        lc.x[i] <- x2[i] + 2/3 * r * cos(phi)
        lc.y[i] <- y2[i] + 2/3 * r * sin(phi)
      }  else {
        spl <- xspline(x = c(c.x1[i], spx[i], c.x2[i]), 
                       y = c(c.y1[i], spy[i], c.y2[i]), shape = 1, 
                       draw = FALSE)
        lines(spl, lwd = sh.lwd[i], col = sh.col[i], 
              lty = sh.lty[i])
        if (code %in% c(2, 3)) {
          x1[i] <- spl$x[3 * length(spl$x)/4]
          y1[i] <- spl$y[3 * length(spl$y)/4]
        }
        if (code %in% c(1, 3)) {
          x2[i] <- spl$x[length(spl$x)/4]
          y2[i] <- spl$y[length(spl$y)/4]
        }
        lc.x[i] <- spl$x[2/3 * length(spl$x)]
        lc.y[i] <- spl$y[2/3 * length(spl$y)]
      }
    }
  }
  if (code %in% c(2, 3)) {
    theta <- atan2((by2 - y1) * uin[2], (bx2 - x1) * uin[1])
    Rep <- rep(length(deg.arr), lx)
    p.x2 <- rep(bx2, Rep)
    p.y2 <- rep(by2, Rep)
    ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
    r.arr <- rep(r.arr, lx)
    x <- (p.x2 + r.arr * cos(ttheta)/uin[1])
    y <- (p.y2 + r.arr * sin(ttheta)/uin[2])
    
    if (open) {
      lines(x[-17], y[-17],
            lwd = h.lwd, col = h.col.bo, lty = h.lty)
    } else { 
      
      # Début
      # x[1]
      # y[1]
      # Fin
      # x[length(x)-1]
      # y[length(x)-1]
      
      # Pointe
      xp <- x[1+(length(x)-1)/2]
      yp <- y[1+(length(x)-1)/2]
      # Point de convergence sur la ligne
      xl <- x[1]+(x[length(x)-1]-x[1])/2
      yl <- y[1]+(y[length(x)-1]-y[1])/2
      
      a <- (yl-yp) / (xl-xp)
      b <- yl - a*xl
      
      d <- sqrt((yl-yp)^2+(xl-xp)^2)
      dp <- beautiful.arrow*d
      
      A <- (1+a^2)
      B <- (-2*xp-2*a*yp+2*a*b)
      C <- (xp^2+yp^2-2*b*yp+b^2-dp^2)
      Delta <- B^2-4*A*C
      
      x01 <- (-B-sqrt(Delta))/(2*A)
      x02 <- (-B+sqrt(Delta))/(2*A)
      
      y01 <- a*x01+b
      y02 <- a*x02+b
      
      d01 <- mean(sqrt((x01-x[1])^2+(y01-y[1])^2), sqrt((x01-x[length(x)-1])^2+(y01-y[length(x)-1])^2))
      d02 <- mean(sqrt((x02-x[1])^2+(y02-y[1])^2), sqrt((x02-x[length(x)-1])^2+(y02-y[length(x)-1])^2))
      
      if (d01<d02) {
        x0 <- x01
        y0 <- y01
      } else {
        x0 <- x02
        y0 <- y02
      }
      
      polygon(c(x0, x[-17], x0), c(y0, y[-17], y0), 
                     col = h.col, lwd = h.lwd, 
                     border = h.col.bo, lty = h.lty)
    }
  }
  if (code %in% c(1, 3)) {
    x1 <- bx1
    y1 <- by1
    tmp <- x1
    x1 <- x2
    x2 <- tmp
    tmp <- y1
    y1 <- y2
    y2 <- tmp
    theta <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
    lx <- length(x1)
    Rep <- rep(length(deg.arr), lx)
    p.x2 <- rep(x2, Rep)
    p.y2 <- rep(y2, Rep)
    ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
    r.arr <- rep(r.arr, lx)
    x <- (p.x2 + r.arr * cos(ttheta)/uin[1])
    y <- (p.y2 + r.arr * sin(ttheta)/uin[2])
    
    if (open) {
      lines(x[-17], y[-17],
            lwd = h.lwd, col = h.col.bo, lty = h.lty)
    } else {
      # Début
      # x[1]
      # y[1]
      # Fin
      # x[length(x)-1]
      # y[length(x)-1]
      
      # Pointe
      xp <- x[1+(length(x)-1)/2]
      yp <- y[1+(length(x)-1)/2]
      # Point de convergence sur la ligne
      xl <- x[1]+(x[length(x)-1]-x[1])/2
      yl <- y[1]+(y[length(x)-1]-y[1])/2
      
      a <- (yl-yp) / (xl-xp)
      b <- yl - a*xl
      
      d <- sqrt((yl-yp)^2+(xl-xp)^2)
      dp <- beautiful.arrow*d
      
      A <- (1+a^2)
      B <- (-2*xp-2*a*yp+2*a*b)
      C <- (xp^2+yp^2-2*b*yp+b^2-dp^2)
      Delta <- B^2-4*A*C
      
      x01 <- (-B-sqrt(Delta))/(2*A)
      x02 <- (-B+sqrt(Delta))/(2*A)
      
      y01 <- a*x01+b
      y02 <- a*x02+b
      
      d01 <- mean(sqrt((x01-x[1])^2+(y01-y[1])^2), sqrt((x01-x[length(x)-1])^2+(y01-y[length(x)-1])^2))
      d02 <- mean(sqrt((x02-x[1])^2+(y02-y[1])^2), sqrt((x02-x[length(x)-1])^2+(y02-y[length(x)-1])^2))
      
      if (d01<d02) {
        x0 <- x01
        y0 <- y01
      } else {
        x0 <- x02
        y0 <- y02
      }
      
      polygon(c(x0, x[-17], x0), c(y0, y[-17], y0), 
              col = h.col, lwd = h.lwd, 
              border = h.col.bo, lty = h.lty)
    }
  }
  invisible(list(lab.x = lc.x, lab.y = lc.y))
}

