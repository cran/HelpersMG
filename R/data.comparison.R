#' data.comparison compares series of data using Akaike weight.
#' @title Data series comparison using Akaike weight
#' @author Marc Girondot
#' @return The probability that a single proportion model is sufficient to explain the data
#' @param ... Series of data (at least two or data are in a table with series in different rows)
#' @param criterion Which criterion is used for model selection. can be AIC, AICc or BIC
#' @param var.equal Should the variances of all series being equal? Default TRUE
#' @description This function is used as a replacement of t.test() to not use p-value. 
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' A <- rnorm(100, 10, 2)
#' B <- rnorm(100, 11.1, 2)
#' data.comparison(A, B, criterion = "BIC", var.equal=TRUE)
#' B <- B[1:10]
#' data.comparison(A, B, criterion = "BIC", var.equal=TRUE)
#' A <- rnorm(100, 10, 2)
#' B <- rnorm(100, 10.1, 2)
#' C <- rnorm(100, 10.5, 2)
#' data.comparison(A, B, criterion = "BIC", var.equal=TRUE)
#' B <- B[1:10]
#' data.comparison(A, B, criterion = "BIC", var.equal=TRUE)
#' t.test(A, B, var.equal=TRUE)
#' 
#' ###########################
#' 
#' replicates <- 500
#' x <- seq(from=8, to=13, by=0.1)
#' dataf <- matrix(rep(NA, replicates*length(x)), nrow=replicates)
#' dataf2 <- dataf
#' 
#' for (repl in 1:replicates) {
#'   
#'   pv <- NULL
#'   aw <- NULL
#'   
#'   for (meanB in x) {
#'     A <- rnorm(100, mean=10, sd=2)
#'     B <- rnorm(100, mean=meanB, sd=4)
#'     
#'     pv <- c(pv, t.test(A, B, var.equal = FALSE)$p.value)
#'     aw <- c(aw, data.comparison(A, B, criterion="BIC", var.equal = FALSE))
#'   }
#'   
#'   dataf[repl, ] <- pv
#'   dataf2[repl, ] <- aw
#'   
#' }
#' 
#' y <- colSums(dataf)/replicates
#' plot(x=x, y=y, type="l", 
#'      bty="n", las=1, xlab="Mean B value (SD = 4)", ylab="Probability", ylim=c(0,1), 
#'      main="Comparison of two data series, A and B (n=100)")
#' y2 <- colSums(dataf2)/replicates
#' lines(x=x, y=y2, type="l", col="red")
#' 
#' x[which(y2>0.05)[1]]
#' rev(x)[which(rev(y2)>0.05)[1]]
#' 
#' x[which(y>0.05)[1]]
#' rev(x)[which(rev(y)>0.05)[1]]
#' 
#' segments(x0=10-1.96*2/10, x1=10+1.96*2/10, y0=1, y1=1, lwd=2)
#' segments(x0=10, x1=10, y0=0.95, y1=1.05, lwd=2)
#' par(xpd=TRUE)
#' text(x=10, y=1.1, labels = "Mean A = 10, SD = 2")
#' 
#' legend("topright", legend=c("p-value\nbased on t-test", 
#'                             "Akaike weight\nfor same A and B model"), 
#'        y.intersp = 2, 
#'        col=c("black", "red"), bty="n", lty=1)
#' 
#' segments(x0=min(x), x1=max(x), y0=0.05, y1=0.05, lty=2)
#' par(xpd = TRUE)
#' text(x=13.05, y=0.05, labels = "0.05", pos=4)
#' 
#' 
#' 
#' ####################################
#' replicates <- 100
#' x <- seq(from=8, to=13, by=0.1)
#' dataf <- matrix(rep(NA, replicates*length(x)), nrow=replicates)
#' dataf2 <- dataf
#' 
#' for (repl in 1:replicates) {
#'   
#'   pv <- NULL
#'   aw <- NULL
#'   
#'   for (meanB in x) {
#'     A <- rnorm(100, 10, 2)
#'     B <- rnorm(100, mean=meanB, sd=2)
#'     
#'     pv <- c(pv, t.test(A, B, var.equal = TRUE)$p.value)
#'     aw <- c(aw, data.comparison(A, B, criterion="BIC", var.equal = TRUE))
#'   }
#'   
#'   dataf[repl, ] <- pv
#'   dataf2[repl, ] <- aw
#'   
#' }
#' 
#' y <- colSums(dataf)/replicates
#' plot(x=x, y=y, type="l", 
#'      bty="n", las=1, xlab="Mean B value (SD = 2)", ylab="Probability", ylim=c(0,1), 
#'      main="Comparison of data series, A and B (n=100)")
#' y2 <- colSums(dataf2)/replicates
#' lines(x=x, y=y2, type="l", col="red")
#' 
#' x[which(y2>0.05)[1]]
#' rev(x)[which(rev(y2)>0.05)[1]]
#' 
#' x[which(y>0.05)[1]]
#' rev(x)[which(rev(y)>0.05)[1]]
#' 
#' 
#' segments(x0=10-1.96*2/10, x1=10+1.96*2/10, y0=1, y1=1, lwd=2)
#' segments(x0=10, x1=10, y0=0.95, y1=1.05, lwd=2)
#' par(xpd=TRUE)
#' text(x=10, y=1.1, labels = "Mean A = 10, SD = 2")
#' 
#' legend("topright", legend=c("p-value\nbased on t-test", 
#'                             "Akaike weight\nfor same A and B model"), 
#'        y.intersp = 2, 
#'        col=c("black", "red"), bty="n", lty=1)
#' 
#' segments(x0=min(x), x1=max(x), y0=0.05, y1=0.05, lty=2)
#' par(xpd=TRUE)
#' text(x=13.05, y=0.05, labels = "0.05", pos=4)
#' }
#' @export


data.comparison <- function(..., criterion = c("BIC", "AIC", "AICc"), var.equal = TRUE) {
  data <- list(...)
  if (length(data) == 1) {
    datax <- data[[1]]
    data <- list()
    for (i in 1:nrow(datax)) data <- c(data, list(datax[i, ]))
  }
  
  if (length(data) == 1) stop("At least two series must be used for comparison")
  
  # data <- list(A, B)
  n <- length(data)
  dtot <- unlist(data)
  sdtot <- sd(dtot)
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
    BIC_s <- BIC_s -2*sum(dnorm(data[[nec]], mean=mean(data[[nec]]), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
    AIC_s <- AIC_s -2*sum(dnorm(data[[nec]], mean=mean(data[[nec]]), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
    AICc_s <- AICc_s -2*sum(dnorm(data[[nec]], mean=mean(data[[nec]]), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
    BIC <- BIC -2*sum(dnorm(data[[nec]], mean=mean(dtot), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
    AIC <- AIC -2*sum(dnorm(data[[nec]], mean=mean(dtot), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
    AICc <-AICc -2*sum(dnorm(data[[nec]], mean=mean(dtot), sd=ifelse(var.equal, sdtot, sd(data[[nec]])), log=TRUE))
  }
  
  
  AIC_s <- AIC_s + 2*nk_AIC_s
  AIC <- AIC + 2*nk_AIC
  BIC_s <- BIC_s + nk_BIC_s*log(length(dtot))
  BIC <- BIC + nk_BIC*log(length(dtot))
  AICc <- AICc + 2*nk_AICc + (2*nk_AICc*(nk_AICc+1))/(length(dtot)-nk_AICc-1)
  AICc_s <- AICc_s + 2*nk_AICc_s + (2*nk_AICc_s*(nk_AICc_s+1))/(length(dtot)-nk_AICc_s-1)
  
  
  crit <- NULL
  AICmin <- min(c(AIC, AIC_s))
  AICcmin <- min(c(AICc, AICc_s))
  BICmin <- min(c(BIC, BIC_s))
  
  if (any(criterion=="AIC")) crit <- c(AW_AIC=exp(-0.5*(AIC-AICmin))/(exp(-0.5*(AIC-AICmin))+exp(-0.5*(AIC_s-AICmin)))         )
  if (any(criterion=="BIC")) crit <- c(crit, AW_BIC=exp(-0.5*(BIC-BICmin))/(exp(-0.5*(BIC-BICmin))+exp(-0.5*(BIC_s-BICmin)))         )
  if (any(criterion=="AICc")) crit <- c(crit, AW_AICc=exp(-0.5*(AICc-AICcmin))/(exp(-0.5*(AICc-AICcmin))+exp(-0.5*(AICc_s-AICcmin)))          )
  
  return(crit)
}

