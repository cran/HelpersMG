#' ExtractAIC.glm returns AIC, AICc or BIC from a glm object
#' @title Return AIC, AICc or BIC from a glm object
#' @author Modified from stats:::extract.AIC.glm
#' @return A numeric named vector of length 2, with first and second elements giving\cr
#' edf	the ‘equivalent degrees of freedom’ for the fitted model fit.\cr
#' x	the Information Criterion for fit.
#' @param fit fitted model, the result of a fitter glm.
#' @param scale unused for glm.
#' @param k numeric specifying the ‘weight’ of the equivalent degrees of freedom (=: edf) part in the AIC formula.
#' @param ... further arguments (currently unused because some functions using this function do not use them).
#' @description For glm fits the family's aic() function is used to compute the AIC.\cr
#' The choice between different criteria is done by setting a global option AIC. It can be checked using show.option=TRUE.
#' Indeed, it is not possible to use the ... parameter due to a bug in some functions of MASS package. If you want to use this function as a replacement for setpAIC(), do extractAIC.glm <- ExtractAIC.glm before.
#' @family AIC
#' @examples
#' extractAIC.glm <- ExtractAIC.glm
#' n <- 100
#' x <- rnorm(n, 20, 2)
#' A <- rnorm(n, 20, 5)
#' g <- glm(x ~ A)
#' extractAIC(g, show.option=TRUE)
#' options(AIC="AIC")
#' extractAIC(g)
#' options(AIC="BIC")
#' extractAIC(g)
#' options(AIC="AICc")
#' extractAIC(g)
#' @export

ExtractAIC.glm <- function (fit, scale = 0, k = 2, ...) {
  p3p <- tryCatch(list(...), error=function(e) list())
  # print(str(p3p))
  type <- getOption("AIC")
  if (is.null(type)) type <- "AIC (default)"
  if (!identical(p3p, list())) {
    if (p3p[["show.option"]]) print(type)
  }
  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  aic <- fit$aic
  if (((type == "AIC (default)")) | (type == "AIC")) return(c(edf, AIC=aic + (k - 2) * edf))
  if (type == "BIC") return(c(edf, BIC=aic + (log(n) - 2) * edf))
  if (type == "AICc") return(c(edf, AICc=aic + (2*edf*(edf+1))/(n-edf-1)))
  if (type == "QAIC") return(c(edf, QAIC=NA))
  if (type == "AICc") return(c(edf, QAICc=NA))
  if (type == "TIC") return(c(edf, TIC=NA))
  return(c(edf, NA))
}

