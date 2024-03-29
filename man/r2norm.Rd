% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2norm.R
\name{r2norm}
\alias{r2norm}
\title{Random generation for Gaussian distributions different at left and right}
\usage{
r2norm(n, mean = 0, sd_low = 1, sd_high = 1)
}
\arguments{
\item{n}{number of observations.}

\item{mean}{vector of means}

\item{sd_low}{vector of standard deviations below the mean.}

\item{sd_high}{vector of standard deviations above the mean.}
}
\value{
r2norm returns random numbers
}
\description{
Random generation for Gaussian distributions different at left and right
}
\details{
r2norm returns random numbers for Gaussian distributions different at left and right
}
\examples{
\dontrun{
n <- r2norm(1000, mean=25, sd_low=2, sd_high=10)

hist(n)
}
}
\seealso{
Other Distributions: 
\code{\link{cutter}()},
\code{\link{dSnbinom}()},
\code{\link{dbeta_new}()},
\code{\link{dcutter}()},
\code{\link{dggamma}()},
\code{\link{logLik.cutter}()},
\code{\link{plot.cutter}()},
\code{\link{print.cutter}()},
\code{\link{rcutter}()},
\code{\link{rmnorm}()},
\code{\link{rnbinom_new}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{Distributions}
