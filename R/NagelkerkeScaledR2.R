#' NagelkerkeScaledR2 returns the scaled R2 defined by Nagelkerke (1991)
#' @title Return the scaled R2 defined by Nagelkerke (1991)
#' @author Marc Girondot
#' @return The scaled R2 value
#' @param x The number of observations
#' @param size Number of trials
#' @param prediction Prediction of x/size
#' @param scaled If TRUE, return the scaled R2
#' @description Return the scaled R2 of a binomial model based on:\cr
#' Nagelkerke NJD (1991) A note on a general definition of the coefficient 
#' of determination. Biometrika 78:691-192.\cr
#' This definition of scaled R2 by Nagelkerke (1991) has the following properties:\cr
#' (i) It is consistent with classical R2, that is the general definition applied to e.g. linear regression yields the classical R2.\cr
#' (ii) It is consistent with maximum likelihood as an estimation method, i.e. the maximum likelihood estimates of the model parameters maximize R2.\cr
#' (iii) It is asymptotically independent of the sample size n.\cr
#' (iv) 1-R2 has the interpretation of the proportion of unexplained 'variation'.\cr
#' (v) It is dimensionless, i.e. it does not depend on the units used.\cr
#' The reported value is similar to the value estimated with nagelkerke() function from rcompanion package but not
#' from the NagelkerkeR2() function from fmsb package. I don't know why.
#' @examples
#' x <- c(10, 9, 6, 4, 3, 1, 0)
#' size <- c(10, 10, 10, 10, 10, 10, 10)
#' prediction <- c(0.9, 0.8, 0.7, 0.5, 0.4, 0.3, 0.2)
#' NagelkerkeScaledR2(x, size, prediction)
#' 
#' # Using the example in fmsb::NagelkerkeR2
#' res <- glm(cbind(ncases,ncontrols) ~ agegp+alcgp+tobgp, data=esoph, family=binomial())
#' NagelkerkeScaledR2(x=esoph$ncases, size = esoph$ncases+esoph$ncontrols, 
#'                    prediction = res$fitted.values)
#' @export


NagelkerkeScaledR2 <- function(x, size, prediction, scaled=TRUE) {
  if (!((length(x) == length(size)) & (length(size) == length(prediction)))) stop("Vectors x, size and prediction must be of same size")
  if (any(prediction > 1) | any(prediction < 0)) stop("Predictions for a binomial distribution must be between 0 and 1.")
  L0 <- sum(dbinom(x,size,rep(sum(x)/sum(size), length(prediction)), log=TRUE))
  Lpar <- sum(dbinom(x, size, prediction, log=TRUE))
  n <- length(size)
  r2 <- 1-(exp(L0-Lpar))^(2/n)
  r2max <- 1-(exp(L0))^(2/n)
  if (!scaled) r2max <- 1
  return(r2/r2max)
}
