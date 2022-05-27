#' RectangleRegression performs rectangle regression
#' @title Return parameters of rectangle regression
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with parameters of rectangle regression
#' @param x1 The first series of data
#' @param x2 The second series of data
#' @param replicate Number of replicates for bootstrap
#' @param x1new Values for x1 to generate x2
#' @description Fit a line using least rectangle method.
#' @family Rectangle Regression
#' @examples
#' x1 <- runif(100, min=10, max=20)
#' x2 <- runif(100, min=10, max=20)+x1
#' 
#' rectreg <- RectangleRegression(x1, x2)
#' 
#' plot(x=x1, y=x2, bty="n", las=1, xlim=c(10, 20), ylim=c(20, 40))
#' abline(a=rectreg$par["Intercept"], b=rectreg$par["Slope"], lwd=2)
#' par(xpd=FALSE)
#' lines(rectreg$x2new["x1new", ], rectreg$x2new["50%", ])
#' lines(rectreg$x2new["x1new", ], rectreg$x2new["2.5%", ], lty=2)
#' lines(rectreg$x2new["x1new", ], rectreg$x2new["97.5%", ], lty=2)
#' 
#' abline(a=rectreg$Intercept[1], b=rectreg$Slope[3], col="red")
#' abline(a=rectreg$Intercept[3], b=rectreg$Slope[1], col="red")
#' 
#' @export


RectangleRegression <- function(x1, x2, replicate=1000, 
                                x1new=seq(from=min(x1), to=max(x1), length.out = 100)) {
  
  frect <- function(x,y){ # droite des moindres rectangles 
    b <- sqrt(var(y)/var(x))*sign(cor(x,y)) 
    r <- c(b, mean(y)-b*mean(x)) 
    names(r) <- c("Slope", "Intercept")
    return(r)
  }

  
  m <- matrix(data=NA, nrow=replicate, ncol=length(x1new))
  c <- NULL
  d <- NULL
  rrx <- frect(x1, x2)
  for (r in 1:replicate) {
    s <- sample(1:length(x1), length(x1), replace = TRUE)
    rr <- frect(x1[s], x2[s])
    c <- c(c, rr[1])
    d <- c(d, rr[2])
    # Prediction
    x2new <- cbind(x1new, 1) %*% rr
    m[r, ] <- x2new
  }
  q <- apply(m, MARGIN=2, FUN=quantile, probs=c(0.025, 0.5, 0.975))
  q <- rbind(q, x1new=x1new)
  qb <- quantile(c, probs=c(0.025, 0.5, 0.975))
  qa <- quantile(d, probs=c(0.025, 0.5, 0.975))
  
  return(list(par=rrx,
              x2new=q,
              Intercept=qa,
              Slope=qb))
}


