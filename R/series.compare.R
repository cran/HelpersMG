#' series.compare compares series of data using Akaike weight.
#' @title Data series comparison using Akaike weight
#' @author Marc Girondot
#' @return The probability that a single proportion model is sufficient to explain the data
#' @param ... Series of data (at least two or data are in a table with series in different rows)
#' @param criterion Which criterion is used for model selection. can be AIC, AICc or BIC
#' @param var.equal Should the variances of all series being equal? Default TRUE
#' @description This function is used as a replacement of t.test() to not use p-value. 
#' @references Girondot, M., Guillon, J.-M., 2018. The w-value: An alternative to t- and X2 tests. Journal of Biostatistics & Biometrics 1, 1-4.
#' @family w-value functions
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' A <- rnorm(100, 10, 2)
#' B <- rnorm(100, 11.1, 2)
#' series.compare(A, B, criterion = "BIC", var.equal=TRUE)
#' B <- B[1:10]
#' series.compare(A, B, criterion = "BIC", var.equal=TRUE)
#' A <- rnorm(100, 10, 2)
#' B <- rnorm(100, 10.1, 2)
#' C <- rnorm(100, 10.5, 2)
#' series.compare(A, B, C, criterion = "BIC", var.equal=TRUE)
#' B <- B[1:10]
#' series.compare(A, B, criterion = "BIC", var.equal=TRUE)
#' t.test(A, B, var.equal=TRUE)
#' # Example with a data.frame
#' series.compare(t(data.frame(A=c(10, 27, 19, 20, NA), B=c(10, 20, NA, NA, NA))))
#' # Test in the context of big data
#' A <- rnorm(10000, 10, 2)
#' B <- rnorm(10000, 10.1, 2)
#' series.compare(A, B, criterion = "BIC", var.equal=TRUE)
#' t.test(A, B, var.equal=TRUE)
#' ###########################
#' w <- NULL
#' p <- NULL
#' 
#' for (i in 1:1000) {
#'   
#'   A <- rnorm(50000, 10, 2)
#'   B <- rnorm(50000, 10.01, 2)
#'   w <- c(w, unname(series.compare(A, B, criterion = "BIC", var.equal=TRUE)[1]))
#'   p <- c(p, t.test(A, B, var.equal=TRUE)$p.value)
#' 
#' }
#' 
#' layout(mat = 1:2)
#' par(mar=c(4, 4, 1, 1)+0.4)
#' hist(p, main="", xlim=c(0, 1), las=1, breaks = (0:20)/20, 
#'      freq=FALSE, xlab = expression(italic("p")*"-value"))
#' hist(w, main="", xlim=c(0, 1), las=1, breaks = (0:20)/20, 
#'     freq=FALSE, xlab = expression(italic("w")*"-value"))
#' ###########################
#' 
#' x <- seq(from=8, to=13, by=0.1)
#' 
#' pv <- NULL
#' aw <- NULL
#' A <- rnorm(100, mean=10, sd=2)
#' B <- A-2
#' 
#' for (meanB in x) {
#'   pv <- c(pv, t.test(A, B, var.equal = FALSE)$p.value)
#'   aw <- c(aw, series.compare(A, B, criterion="BIC", var.equal = FALSE)[1])
#'   B <- B + 0.1
#' }
#' 
#' par(mar=c(4, 4, 2, 1)+0.4)
#' y <- pv
#' plot(x=x, y=y, type="l", lwd=2,
#'      bty="n", las=1, xlab="Mean B value (SD = 4)", ylab="Probability", ylim=c(0,1), 
#'      main="")
#' y2 <- aw
#' lines(x=x, y=y2, type="l", col="red", lwd=2)
#' 
# w-value
#' l1 <- which(aw>0.05)[1]
#' l2 <- max(which(aw>0.05))
#' 
#' aw[l1]
#' pv[l1]
#' 
#' aw[l2]
#' pv[l2]
#' 
# p-value
#' l1 <- which(pv>0.05)[1]
#' l2 <- max(which(pv>0.05))
#' 
#' aw[l1]
#' pv[l1]
#' 
#' aw[l2]
#' pv[l2]
#' 
#' par(xpd=TRUE)
#' segments(x0=10-1.96*2/10, x1=10+1.96*2/10, y0=1.1, y1=1.1, lwd=2)
#' segments(x0=10, x1=10, y0=1.15, y1=1.05, lwd=2)
#' par(xpd=TRUE)
#' text(x=10.5, y=1.1, labels = "Mean A = 10, SD = 2", pos=4)
#' 
#' v1 <- c(expression(italic("p")*"-value"), expression("based on "*italic("t")*"-test"))
#' v2 <- c(expression(italic("w")*"-value for A"), expression("and B identical models"))
#' legend("topright", legend=c(v1, v2), 
#'        y.intersp = 1, 
#'        col=c("black", "black", "red", "red"), bty="n", lty=c(1, 0, 1, 0))
#' 
#' segments(x0=min(x), x1=max(x), y0=0.05, y1=0.05, lty=2)
#' par(xpd = TRUE)
#' text(x=13.05, y=0.05, labels = "0.05", pos=4)
#' }
#' @export


series.compare <- function(..., criterion = c("BIC", "AIC", "AICc"), var.equal = TRUE) {
  
  # data <- list(A, B); criterion = "BIC"; var.equal=TRUE
  
  data <- list(...)
  if (length(data) == 1) {
    datax <- data[[1]]
    data <- list()
    for (i in 1:nrow(datax)) data <- c(data, list(datax[i, !is.na(datax[i, ])]))
  }
  
  if (length(data) == 1) stop("At least two series must be used for comparison")
  
  # data <- list(A, B)
  n <- length(data)
  dtot <- unlist(data)
  sdtot <- sd(dtot)
  n_tot <- length(dtot)
  
  AIC_s <- 0
  BIC_s <- 0
  AICc_s <- 0
  BIC <- 0
  AIC <- 0
  AICc <-  0
  
  nk_AIC <- 1+ifelse(var.equal, 1, n)
  nk_AIC_s <- n+ifelse(var.equal, 1, n)
  nk_BIC <- 1+ifelse(var.equal, 1, n)
  nk_BIC_s <- n+ifelse(var.equal, 1, n)
  nk_AICc <- 1+ifelse(var.equal, 1, n)
  nk_AICc_s <- n+ifelse(var.equal, 1, n)
  
  
  for (nec in 1:n) {
    logL <- sum(dnorm(data[[nec]], mean=mean(data[[nec]]), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
    BIC_s <- BIC_s -2*logL
    AIC_s <- AIC_s -2*logL
    AICc_s <- AICc_s -2*logL
    logL <- sum(dnorm(data[[nec]], mean=mean(dtot), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
    BIC <- BIC -2*logL
    AIC <- AIC -2*logL
    AICc <-AICc -2*logL
  }
  
  AIC_s <- AIC_s + 2*nk_AIC_s
  AIC <- AIC + 2*nk_AIC
  BIC_s <- BIC_s + nk_BIC_s*log(n_tot)
  BIC <- BIC + nk_BIC*log(n_tot)
  AICc <- AICc + 2*nk_AICc + (2*nk_AICc*(nk_AICc+1))/(n_tot-nk_AICc-1)
  AICc_s <- AICc_s + 2*nk_AICc_s + (2*nk_AICc_s*(nk_AICc_s+1))/(n_tot-nk_AICc_s-1)
  
  
  crit <- NULL
  AICmin <- min(c(AIC, AIC_s))
  AICcmin <- min(c(AICc, AICc_s))
  BICmin <- min(c(BIC, BIC_s))
  
  if (any(criterion=="AIC")) crit <- c("AICw(identical)"=exp(-0.5*(AIC-AICmin))/(exp(-0.5*(AIC-AICmin))+exp(-0.5*(AIC_s-AICmin))), 
                                       "AICw(different)"=1-(exp(-0.5*(AIC-AICmin))/(exp(-0.5*(AIC-AICmin))+exp(-0.5*(AIC_s-AICmin))))
                                       )
  if (any(criterion=="BIC")) crit <- c(crit, "BICw(identical)"=exp(-0.5*(BIC-BICmin))/(exp(-0.5*(BIC-BICmin))+exp(-0.5*(BIC_s-BICmin))), 
                                                    "BICw(different)"=1-(exp(-0.5*(BIC-BICmin))/(exp(-0.5*(BIC-BICmin))+exp(-0.5*(BIC_s-BICmin))))                  
                                       )
  if (any(criterion=="AICc")) crit <- c(crit, "AICcw(identical)"=exp(-0.5*(AICc-AICcmin))/(exp(-0.5*(AICc-AICcmin))+exp(-0.5*(AICc_s-AICcmin))), 
                                        "AICcw(different)"=1-(exp(-0.5*(AICc-AICcmin))/(exp(-0.5*(AICc-AICcmin))+exp(-0.5*(AICc_s-AICcmin))))
                                        )
  
  return(crit)
}

