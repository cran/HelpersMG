#' invlogit returns the inverse logit
#' @title Return the inverse logit
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A value
#' @param n The value to inverse to get the probability
#' @description Return the inverse logit.
#' @family logit
#' @examples
#' n <- logit(0.5)
#' invlogit(n)
#' @export

invlogit <- function(n) {1/(1+exp(-n))}
