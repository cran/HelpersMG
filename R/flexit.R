#' Return the flexit value
#' @title Return the flexit
#' @author Marc Girondot
#' @return A vector with the probabilities
#' @param x The values at which the flexit model must be calculated
#' @param par The vector with P, S, K1, and K2 values
#' @param P P value
#' @param S S value
#' @param K1 K1 value
#' @param K2 K2 value
#' @param zero Value to replace zero
#' @param error0 Value to return if an error is observed toward 0
#' @param error1 Value to return if an error is observed toward 1
#' @description Return a vector with the probabilities.
#' The flexit equation is not still published :
#' \deqn{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}
#' \deqn{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}
#' with:\cr
#'      \deqn{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}
#'      \deqn{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}
#' @family logit
#' @examples
#' n <- flexit(x=1:100, par=c(P=50, S=0.001, K1=0.01, K2=0.02))
#' n <- flexit(x=1:100, P=50, S=0.001, K1=0.01, K2=0.02)
#' @export

flexit <- function(x, par=NULL, P=NULL, S=NULL, K1=NULL, K2=NULL, zero=1E-9, error0 = 0, error1 = 1) {

  # Peut-être  encore des problèmes de exp(K1 ou K2)
  if (!is.null(par)) {
    K1 <- ifelse(par["K1"] == 0, zero, par["K1"])
    K2 <- ifelse(par["K2"] == 0, zero, par["K2"])
    S <- par["S"]
    P <- par["P"]
  }
  
  if (is.null(K1)) K1 <- zero
  if (is.null(K2)) K2 <- zero
  K1 <- ifelse(K1 == 0, zero, K1)
  K2 <- ifelse(K2 == 0, zero, K2)
  
  if (is.infinite(2^(K1))) {
    S1 <- K1*S
    K1 <- sign(K1)*500
  } else {
    S1 <- (2^(K1 - 1)*K1*S)/(2^(K1) - 1)
  }
  if (is.infinite(2^(K2))) {
    S2 <- K2*S
    K2 <- sign(K2)*500
  } else {
    S2 <- (2^(K2 - 1)*K2*S)/(2^(K2) - 1)
  }
  
  p <- ifelse(x < P,  ifelse((1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))<0,
                                 ifelse(S1 < 0, error1, error0),
                                 (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1))
              ,
              ifelse((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))<0,
                     ifelse(S2 < 0, error0, error1),
                     1-(1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2))
  )
  
  return(p)
}
