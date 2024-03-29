% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setPriors.R
\name{setPriors}
\alias{setPriors}
\title{Set priors for MHalgoGen()}
\usage{
setPriors(
  par = stop("A vector with init values is necessary."),
  se = NULL,
  density = "dunif",
  rules = NULL,
  silent = FALSE
)
}
\arguments{
\item{par}{Named vector with init value of parameters}

\item{se}{Named vector with standard error of parameters}

\item{density}{Named vector with density or single value}

\item{rules}{List of rules to define priors}

\item{silent}{If TRUE, do not show warning.}
}
\value{
Return a data.frame with priors
}
\description{
Set priors for MHalgoGen()\cr
}
\details{
setPriors is a general function to set priors for MHalgoGen()
}
\examples{
\dontrun{
library(HelpersMG)
rules <- rbind(data.frame(Name="^a", Min=0, Max="x*2"), 
              data.frame(Name="^b", Min=0, Max=100))
par <- c(a0=10, a1=2, b2=20)
(p <- setPriors(par=par, se=NULL, density="dgamma", rules=rules))
(p <- setPriors(par=par, se=NULL, density="dnorm", rules=rules))
(p <- setPriors(par=par, se=NULL, density="dunif", rules=rules))
par <- c(a0=10, a1=2, b2=20, b1=-1)
(p <- setPriors(par=par, se=NULL, density="dgamma", rules=rules))
}
}
\seealso{
Other mcmcComposite functions: 
\code{\link{MHalgoGen}()},
\code{\link{as.mcmc.mcmcComposite}()},
\code{\link{as.parameters}()},
\code{\link{as.quantiles}()},
\code{\link{merge.mcmcComposite}()},
\code{\link{plot.PriorsmcmcComposite}()},
\code{\link{plot.mcmcComposite}()},
\code{\link{summary.mcmcComposite}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{mcmcComposite functions}
