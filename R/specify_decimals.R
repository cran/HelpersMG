#' specify_decimals format a number with specified number of decimals
#' @title Return a number as character with specified number of decimals
#' @author Marc Girondot
#' @return A character
#' @param x The numbers to be formated
#' @param decimals Number of decimals to print
#' @param decimal.point Character to be used as decimal point
#' @description Return a number as character with specified number of decimals.
#' @examples
#' specify_decimal(x=pi, decimal.point=".")
#' specify_decimal(x=pi, decimals=4, decimal.point=".")
#' specify_decimal(x=c(pi, exp(1)), decimals=3, decimal.point=",")
#' specify_decimal(x=c(pi, exp(1)), decimal.point=",")
#' specify_decimal(x=c(pi*10, pi, pi/10, pi/100, pi/1000))
#' @export

specify_decimal <- function(x, 
                            decimals=NULL, 
                            decimal.point=".") {
  if (is.character(x)) x <- as.numeric(x)
  if (is.null(decimals)) {
    na <- is.na(x)
    x[is.na(x)] <- 1
    og <- log10(abs(x))
    og <- ifelse(og == -Inf, 1, og)
    decimals <- ifelse(og > 0, 3, abs(floor(og))+3)
  } else {
    na <- rep(FALSE, length(x))
  }
  decimals <- rep(decimals, length(x))[1:length(x)]
  cx <- sapply(seq_along(x), FUN=function(i) formatC(x[i], digits = decimals[i], format = "f"))
  if (any(na)) cx[na] <- "NA"
  return(gsub("\\.", decimal.point, unname(cx)))
}

