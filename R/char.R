#' char returns the characters defined by the codes
#' @title Return the characters defined by the codes
#' @author Based on this blog: http://datadebrief.blogspot.com/2011/03/ascii-code-table-in-r.html
#' @return A string with characters defined by the codes
#' @param n The code to be used to return a character
#' @description Return a string with characters defined by the codes.
#' @family Characters
#' @examples
#' char(65:75)
#' char(unlist(tapply(144:175, 144:175, function(x) {c(208, x)})))
#' @export


char <- function(n) {rawToChar(as.raw(n))}

