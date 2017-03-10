#' contingencyTable.compare compares contingency table using Akaike weight.
#' @title Contingency table comparison using Akaike weight
#' @author Marc Girondot
#' @return The probability that a single proportion model is sufficient to explain the data
#' @param table A matrix or a data.frame with series in rows and number of each category in column
#' @param criterion Which criterion is used for model selection
#' @param probs Series of probabilities used for conformity comparison
#' @description This function is used as a replacement of chisq.test() to not use p-value. 
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' 
#' # Symmetry of Lepidochelys olivacea scutes
#' table <- t(data.frame(SriLanka=c(200, 157), AfricaAtl=c(19, 12), 
#'                       Guyana=c(8, 6), Suriname=c(162, 88), 
#'                       MexicoPac1984=c(42, 34), MexicoPac2014Dead=c(8, 9),
#'                       MexicoPac2014Alive=c(13, 12), 
#'                       row.names =c("Symmetric", "Asymmetric")))
#' table
#' contingencyTable.compare(table)
#' 
#' table <- t(data.frame(SriLanka=c(200, 157), AfricaAtl=c(19, 12), Guyana=c(8, 6),
#'                       Suriname=c(162, 88), MexicoPac1984=c(42, 34), 
#'                       MexicoPac2014Dead=c(8, 9),
#'                       MexicoPac2014Alive=c(13, 12), Lepidochelys.kempii=c(99, 1), 
#'                       row.names =c("Symmetric", "Asymmetric")))
#' table
#' contingencyTable.compare(table)
#' 
#' # Conformity to a model
#' table <- matrix(c(33, 12, 25, 75), ncol = 2, byrow = TRUE)
#' probs <- c(0.5, 0.5)
#' contingencyTable.compare(table, probs=probs)
#' 
#' # Conformity to a model
#' table <- matrix(c(33, 12), ncol = 2, byrow = TRUE)
#' probs <- c(0.5, 0.5)
#' contingencyTable.compare(table, probs=probs)
#' 
#' # Conformity to a model
#' table <- matrix(c(33, 12, 8, 25, 75, 9), ncol = 3, byrow = TRUE)
#' probs <- c(0.8, 0.1, 0.1)
#' contingencyTable.compare(table, probs=probs)
#' 
#' # Comparison of chisq.test() and this function
#' table <- matrix(c(NA, NA, 25, 75), ncol = 2, byrow = TRUE)
#' 
#' pv <- NULL
#' aw <- NULL
#' par(new=FALSE)
#' n <- 100
#' 
#' for (GroupA in 0:n) {
#'   table[1, 1] <- GroupA
#'   table[1, 2] <- n-GroupA
#'   pv <- c(pv, chisq.test(table)$p.value)
#'   aw <- c(aw, contingencyTable.compare(table, criterion="BIC")[1])
#' }
#' 
#' x <- 0:n
#' y <- pv
#' y2 <- aw
#' plot(x=x, y=y, type="l", bty="n", las=1, xlab="Number of type P in Group B", ylab="Probability", 
#'      main="", lwd=2)
#' lines(x=x, y=y2, type="l", col="red", lwd=2)
#' 
#' # w-value
#' (l1 <- x[which(aw>0.05)[1]])
#' (l2 <- rev(x)[which(rev(aw)>0.05)[1]])
#' 
#' aw[l1]
#' pv[l1]
#' 
#' aw[l2+2]
#' pv[l2+2]
#' 
#' # p-value
#' l1 <- which(pv>0.05)[1]
#' l2 <- max(which(pv>0.05))
#' 
#' aw[l1]
#' pv[l1]
#' 
#' aw[l2]
#' pv[l2]
#' 
#' y[which(y2>0.05)[1]]
#' y[which(rev(y2)>0.05)[1]]
#' 
#' par(xpd=TRUE)
#' text(x=25, y=1.15, labels="Group A: 25 type P / 100", pos=1)
#' 
#' segments(x0=25, y0=0, x1=25, y1=1, lty=3)
#' 
#' # plot(1, 1)
#' 
#' v1 <- c(expression(italic("p")*"-value"), expression("after "*chi^2*"-test"))
#' v2 <- c(expression(italic("w")*"-value for A"), expression("and B identical models"))
#' legend("topright", legend=c(v1, v2), 
#'        y.intersp = 1, 
#'        col=c("black", "black", "red", "red"), bty="n", lty=c(1, 0, 1, 0))
#' 
#' segments(x0=0, x1=n, y0=0.05, y1=0.05, lty=2)
#' text(x=101, y=0.05, labels = "0.05", pos=4)
#' }
#' @export


contingencyTable.compare <- function(table, criterion=c("AIC", "AICc", "BIC"), probs=NULL) {
  # table <- matrix(c(10, 12, 19, 21, 7, 8), ncol = 2)
  
  
  nseries <- nrow(table)
  c <- ncol(table)
  k <- c-1
  n <- nseries*k
  
  AIC_s <- 0
  AIC <- 0
  BIC_s <- 0
  BIC <- 0
  AICc_s <- 0
  AICc <- 0
  for (nec in 1:nseries){
    l <- table[nec, ]
    if (is.null(probs)) {
      prob <- l/sum(l)
    } else {
      prob <- probs
    }
    logL <- dmultinom(l, size = sum(l), prob = prob, log = TRUE)
    BIC_s <- BIC_s -2*logL
    AIC_s <- AIC_s -2*logL
    AICc_s <- AICc_s -2*logL
    if (is.null(probs)) {
      prob <- colSums(table)/sum(colSums(table))
    } else {
      prob <- probs
    }
    logL <- dmultinom(l, size = sum(l), prob = prob, log = TRUE)
    AIC <- AIC -2*logL
    BIC <- BIC -2*logL
    AICc <- AICc -2*logL
  }
  
  k_common <- ifelse(is.null(probs), k, 0)
  k_s <- n
    
  AIC <- AIC + 2*k_common
  AIC_s <- AIC_s + 2*k_s
  
  BIC <- BIC + k_common*log(n)
  BIC_s <- BIC_s + k_s*log(n)
  
  AICc <- AICc + 2*k_common + (2*k_common*(k_common+1))/(n-k_common-1)
  AICc_s <- AICc_s + 2*k_s + (2*k_s*(k_s+1))/(n-k_s-1)
  
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
  
  if (any(criterion=="AICc")) {
   if ((2*k_s*(k_s+1))/(n-k_s-1)<0) {
      warning("The AICc cannot be estimated.")
      crit <- c(crit, "AICcw(identical)"=NA, "AICcw(different)"=NA)
    } else {
      crit <- c(crit, "AICcw(identical)"=exp(-0.5*(AICc-AICcmin))/(exp(-0.5*(AICc-AICcmin))+exp(-0.5*(AICc_s-AICcmin))), 
                                        "AICcw(different)"=1-(exp(-0.5*(AICc-AICcmin))/(exp(-0.5*(AICc-AICcmin))+exp(-0.5*(AICc_s-AICcmin))))
                )
    }
  }
    
  
  return(crit)
}

