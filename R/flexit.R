#' Return the flexit value
#' @title Return the flexit
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A vector with the probabilities
#' @param x The values at which the flexit model must be calculated
#' @param par The vector with P, S, K1, and K2 values
#' @param P P value
#' @param S S value
#' @param K1 K1 value
#' @param K2 K2 value
#' @param Min Min value for scaled flexit model
#' @param Max Max value for scaled flexit model
#' @param zero Value to replace zero
#' @param error0 Value to return if an error is observed toward 0
#' @param error1 Value to return if an error is observed toward 1
#' @description Return a vector with the probabilities.
#' The flexit equation is published in:\cr
#' Abreu-Grobois, F.A., Morales-Mérida, B.A., Hart, C.E., Guillon, J.-M., Godfrey, M.H., 
#' Navarro, E. & Girondot, M. (2020) Recent advances on the estimation of the thermal 
#' reaction norm for sex ratios. PeerJ, 8, e8451.\cr
#' If dose < P then \eqn{(1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}{(1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}\cr
#' If dose > P then \eqn{1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}{1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}\cr
#' with:\cr
#'      \deqn{S1 = (2^(K1 - 1) * S * K1)/(2^K1 - 1)}{S1 = (2^(K1 - 1) * S * K1)/(2^K1 - 1)}
#'      \deqn{S2 = (2^(K2 - 1) * S * K2)/(2^K2 - 1)}{S2 = (2^(K2 - 1) * S * K2)/(2^K2 - 1)}
#' \cr
#' If \eqn{2^K1}{2^K1} is too large to be estimated, the approximation \eqn{S1 = S*K1/2}{S1 = S*K1/2} is used.\cr
#' Demonstration:\cr
#' \deqn{S1 = (2^(K1 - 1) * S * K1)/(2^K1 - 1)}{S1 = (2^(K1 - 1) * S * K1)/(2^K1 - 1)}
#' \deqn{S1 = exp(log((2^(K1 - 1) * S * K1)/(2^K1 - 1)))}{S1 = exp(log((2^(K1 - 1) * S * K1)/(2^K1 - 1)))}
#' \deqn{S1 = exp(log(2^(K1 - 1)) + log(S * K1) - log(2^K1 - 1))}{S1 = exp(log(2^(K1 - 1)) + log(S * K1) - log(2^K1 - 1))}
#' When \eqn{K1}{K1} is very large, \eqn{2^K1 - 1 = 2^K1}{2^K1 - 1 = 2^K1} then
#' \deqn{S1 = exp((K1 - 1) * log(2) + log(S * K1) - K1 * log(2))}{S1 = exp((K1 - 1) * log(2) + log(S * K1) - K1 * log(2))}
#' \deqn{S1 = exp((K1 * log(2) - log(2) + log(S * K1) - K1 * log(2))}{S1 = exp((K1 * log(2) - log(2) + log(S * K1) - K1 * log(2))}
#' \deqn{S1 = exp(log(S * K1)- log(2))}{S1 = exp(log(S * K1)- log(2))}
#' \deqn{S1 = S * K1 / 2}{S1 = S * K1 / 2}
#' If \eqn{2^K2}{2^K2} is too large to be estimated, the approximation \eqn{S2 = S*K2/2}{S2 = S*K2/2} is used.\cr
#' If \eqn{(1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}{(1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)} is not finite, 
#' the following approximation is used:\cr
#' \deqn{exp((-1/K1)*(K1*log(2)+(4*S1*(P-x))))}{exp((-1/K1)*(K1*log(2)+(4*S1*(P-x))))}
#' If \eqn{1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}{1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)} is not finite, 
#' the following approximation is used:\cr
#' \deqn{1 - exp((-1/K2)*(K2*log(2)+(4*S2*(x - P))))}{1 - exp((-1/K2)*(K2*log(2)+(4*S2*(x - P))))}
#' @family logit
#' @examples
#' n <- flexit(x=1:100, par=c(P=50, S=0.001, K1=0.01, K2=0.02))
#' n <- flexit(x=1:100, P=50, S=0.001, K1=0.01, K2=0.02)
#' 
#' 1/(1+exp(0.01*4*(50-1:100)))
#' flexit(1:100, P=50, S=0.01, K1=1, K2=1)
#' @export

flexit <- function(x, par=NULL, P=NULL, S=NULL, K1=NULL, K2=NULL, Min=0, Max=1, zero=1E-9, error0 = 0, error1 = 1) {

  # par=NULL; P=NULL; S=NULL; K1=NULL; K2=NULL; Min=0; Max=1; zero=1E-9; error0 = 0; error1 = 1
  
#             par <- c('P' = 0.66325432700364861, 
#               'S' = 6.68963311308306, 
#               'Min' = 0.16643548495805249, 
#               'Max' = 1, 
#               'K1' = 1, 
#               'K2' = 1)
#             x <- seq(from=0.005, to=0.995, by=0.01)
#             zero=1E-9; error0 = 0; error1 = 1
  
  # Peut-être  encore des problèmes de exp(K1 ou K2)
  if (!is.null(par)) {
    if (is.na(par["K1"])) par <- c(par, K1=1)
    if (is.na(par["K2"])) par <- c(par, K2=1)
    if (is.na(par["Min"])) par <- c(par, Min=0)
    if (is.na(par["Max"])) par <- c(par, Max=1)
    K1 <- ifelse(par["K1"] == 0, zero, par["K1"])
    K2 <- ifelse(par["K2"] == 0, zero, par["K2"])
    Min <- par["Min"]
    Max <- par["Max"]
    S <- par["S"]
    P <- par["P"]
  }
  
  lg <- max(c(length(x), length(K1), length(K2), length(P), length(S), length(Min), length(Max)))
  
  if (is.null(Min)) {
    Min <- rep(zero, lg)
  } else {
    if (length(Min) != lg) Min <- rep(Min, lg)[1:lg]
  }
  if (is.null(Max)) {
    Max <- rep(1-zero, lg)
  } else {
    if (length(Max) != lg) Max <- rep(Max, lg)[1:lg]
  }
  if (is.null(K1)) {
    K1 <- rep(1, lg)
  } else {
    if (length(K1) != lg) K1 <- rep(K1, lg)[1:lg]
  }
  if (is.null(K2)) {
    K2 <- rep(1, lg)
  } else {
    if (length(K2) != lg) K2 <- rep(K2, lg)[1:lg]
  }
  K1 <- ifelse(K1 == 0, zero, K1)
  K2 <- ifelse(K2 == 0, zero, K2)
  
  # Je les mets tous à la même taille
  if (length(P) != lg) P <- rep(P, lg)[1:lg]
  if (length(S) != lg) S <- rep(S, lg)[1:lg]
  if (length(x) != lg) x <- rep(x, lg)[1:lg]
  
  S1 <- ifelse(is.infinite(2^(K1)), K1*S/2, (2^(K1 - 1)*K1*S)/(2^(K1) - 1))
  
  # if (is.infinite(2^(K1))) {
  #   S1 <- K1*S/2
  # } else {
  #   S1 <- (2^(K1 - 1)*K1*S)/(2^(K1) - 1)
  # }
  
  Test1 <- (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))
  Test1_p <- ifelse(!is.infinite(Test1), 
                    Test1^(-1/K1), 
                    exp((-1/K1)*(K1*log(2)+(4*S1*(P-x)))))
  
  S2 <- ifelse(is.infinite(2^(K2)), K2*S/2, (2^(K2 - 1)*K2*S)/(2^(K2) - 1))
  
  # if (is.infinite(2^(K2))) {
  #   S2 <- K2*S/2
  # } else {
  #   S2 <- (2^(K2 - 1)*K2*S)/(2^(K2) - 1)
  # }
  
  Test2 <- (1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))
  Test2_p <- ifelse(!is.infinite(Test2), 
                    1 - Test2^(-1/K2), 
                    1 - exp((-1/K2)*(K2*log(2)+(4*S2*(x - P)))))
  
  
  p <- ifelse(x < P,  ifelse(Test1 < 0,
                                 ifelse(S1 < 0, error1, error0),
                             Test1_p)
              ,
              ifelse(Test2 < 0,
                     ifelse(S2 < 0, error0, error1),
                     Test2_p)
  )
  
  p <- (Max - Min) * p + Min
  
  return(unname(p))
}
