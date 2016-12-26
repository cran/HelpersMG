#' table.comparison compares contingency table using Akaike weight.
#' @title Contingency table comparison using Akaike weight
#' @author Marc Girondot
#' @return The probability that a single proportion model is sufficient to explain the data
#' @param table A matrix or a data.frame with series in rows and number of each category in column
#' @param criterion Which criterion is used for model selection. Only AIC is supported
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
#' table.comparison(table)
#' 
#' table <- t(data.frame(SriLanka=c(200, 157), AfricaAtl=c(19, 12), Guyana=c(8, 6),
#'                       Suriname=c(162, 88), MexicoPac1984=c(42, 34), 
#'                       MexicoPac2014Dead=c(8, 9),
#'                       MexicoPac2014Alive=c(13, 12), Lepidochelys.kempii=c(99, 1), 
#'                       row.names =c("Symmetric", "Asymmetric")))
#' table
#' table.comparison(table)
#' 
#' # Comparison of chisq.test() and this function
#' table <- matrix(c(NA, NA, 25, 100), ncol = 2, byrow = TRUE)
#' pv <- NULL
#' aw <- NULL
#' 
#' n <- 100
#' 
#' for (GroupA in 0:n) {
#'   table[1, 1] <- GroupA
#'   table[1, 2] <- n-GroupA
#'   pv <- c(pv, chisq.test(table)$p.value)
#'   aw <- c(aw, table.comparison(table))
#' }
#' 
#' x <- 0:n
#' y <- pv
#' y2 <- aw
#' plot(x=x, y=y, type="l", bty="n", las=1,
#'      xlab="Number in Group A / 100", ylab="Probability", 
#'      main="Comparison of a contingency table with two rows (A and B)")
#' lines(x=x, y=y2, type="l", col="red")
#' 
#' x[which(y2>0.05)[1]]
#' rev(x)[which(rev(y2)>0.05)[1]]
#' 
#' x[which(y>0.05)[1]]
#' rev(x)[which(rev(y)>0.05)[1]]
#' 
#' par(xpd=TRUE)
#' text(x=20, y=1.15, labels="Group B: 25 / 100", pos=1)
#' 
#' legend("topright", legend=c("p-value for chi-2", 
#'                             "Akaike weight for same A and B models"), 
#'        col=c("black", "red"), bty="n", lty=1)
#' segments(x0=0, x1=n, y0=0.05, y1=0.05, lty=2)
#' text(x=101, y=0.05, labels = "0.05", pos=4)
#' }
#' @export


table.comparison <- function(table, criterion="AIC") {
  # table <- matrix(c(10, 12, 19, 21, 7, 8), ncol = 2)
  n <- nrow(table)
  c <- ncol(table)
  k <- c-1
  AIC_s <- 0
  AIC <- 2*k
#  BIC_s <- NULL
#  BIC <- NULL
#  AICc_s <- NULL
#  AICc <- NULL
  for (nec in 1:n){
    l <- table[nec, ]
    logL <- dmultinom(l, size = sum(l), prob = l/sum(l), log = TRUE)
    
    k <- c-1
    # BIC_s <- BIC_s -2*logL + k * log(length(data[[nec]]))
    AIC_s <- AIC_s -2*logL + 2*k
    # AICc_s <- AICc_s -2*logL + 2*k + (2*k*(k+1))/(length(data[[nec]])-k-1)
    logL <- dmultinom(l, size = sum(l), prob = colSums(table)/sum(colSums(table)), log = TRUE)
    k <- 0
    AIC <- AIC -2*logL + 2*k
  }
  crit <- NULL
  AICmin <- min(c(AIC, AIC_s))

  if (any(criterion=="AIC")) crit <- c(AW_AIC=exp(-0.5*(AIC-AICmin))/(exp(-0.5*(AIC-AICmin))+exp(-0.5*(AIC_s-AICmin)))         )
 
  return(crit)
}

