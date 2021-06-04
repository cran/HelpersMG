#' specify_decimals format a number with specified number of decimals
#' @title Return a number as character with specified number of decimals
#' @author Marc Girondot
#' @return A character
#' @param x The numbers to be formated
#' @param decimals Number of decimals to print
#' @param decimal.point Character to be used as decimal point
#' @description Return a number as character with specified number of decimals.
#' @examples
#' specify_decimal(x=pi, decimals=3, decimal.point=".")
#' specify_decimal(x=c(pi, exp(1)), decimals=3, decimal.point=",")
#' @export

specify_decimal <- function(x, 
                            decimals=3, 
                            decimal.point=".") 
  gsub("\\.", decimal.point, unname(formatC(x, digits = decimals, format = "f")))


