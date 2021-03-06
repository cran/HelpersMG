% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexit.R
\name{flexit}
\alias{flexit}
\title{Return the flexit}
\usage{
flexit(
  x,
  par = NULL,
  P = NULL,
  S = NULL,
  K1 = NULL,
  K2 = NULL,
  zero = 1e-09,
  error0 = 0,
  error1 = 1
)
}
\arguments{
\item{x}{The values at which the flexit model must be calculated}

\item{par}{The vector with P, S, K1, and K2 values}

\item{P}{P value}

\item{S}{S value}

\item{K1}{K1 value}

\item{K2}{K2 value}

\item{zero}{Value to replace zero}

\item{error0}{Value to return if an error is observed toward 0}

\item{error1}{Value to return if an error is observed toward 1}
}
\value{
A vector with the probabilities
}
\description{
Return a vector with the probabilities.
The flexit equation is published in:\cr
Abreu-Grobois, F.A., Morales-Mérida, B.A., Hart, C.E., Guillon, J.-M., Godfrey, M.H., 
Navarro, E. & Girondot, M. (2020) Recent advances on the estimation of the thermal 
reaction norm for sex ratios. PeerJ, 8, e8451.\cr
\deqn{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}
\deqn{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}
with:\cr
     \deqn{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}
     \deqn{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}
}
\details{
Return the flexit value
}
\examples{
n <- flexit(x=1:100, par=c(P=50, S=0.001, K1=0.01, K2=0.02))
n <- flexit(x=1:100, P=50, S=0.001, K1=0.01, K2=0.02)
}
\seealso{
Other logit: 
\code{\link{invlogit}()},
\code{\link{logit}()}
}
\author{
Marc Girondot
}
\concept{logit}
