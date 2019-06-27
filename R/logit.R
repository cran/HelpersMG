#' logit returns the logit
#' @title Return the logit
#' @author Marc Girondot
#' @return A value
#' @param p The probability
#' @description Return the logit.
#' @family logit
#' @examples
#' n <- logit(0.5)
#' invlogit(n)
#' @export

logit <- function(p) {log(p/(1-p))}
