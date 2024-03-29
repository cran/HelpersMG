% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NagelkerkeScaledR2.R
\name{NagelkerkeScaledR2}
\alias{NagelkerkeScaledR2}
\title{Return the scaled R2 defined by Nagelkerke (1991)}
\usage{
NagelkerkeScaledR2(x, size, prediction, scaled = TRUE)
}
\arguments{
\item{x}{The number of observations}

\item{size}{Number of trials}

\item{prediction}{Prediction of x/size}

\item{scaled}{If TRUE, return the scaled R2}
}
\value{
The scaled R2 value
}
\description{
Return the scaled R2 of a binomial model based on:\cr
Nagelkerke NJD (1991) A note on a general definition of the coefficient
of determination. Biometrika 78:691-192.\cr
This definition of scaled R2 by Nagelkerke (1991) has the following properties:\cr
(i) It is consistent with classical R2, that is the general definition applied to e.g. linear regression yields the classical R2.\cr
(ii) It is consistent with maximum likelihood as an estimation method, i.e. the maximum likelihood estimates of the model parameters maximize R2.\cr
(iii) It is asymptotically independent of the sample size n.\cr
(iv) 1-R2 has the interpretation of the proportion of unexplained 'variation'.\cr
(v) It is dimensionless, i.e. it does not depend on the units used.\cr
The reported value is similar to the value estimated with nagelkerke() function from rcompanion package but not
from the NagelkerkeR2() function from fmsb package. I don't know why.
}
\details{
NagelkerkeScaledR2 returns the scaled R2 defined by Nagelkerke (1991)
}
\examples{
x <- c(10, 9, 6, 4, 3, 1, 0)
size <- c(10, 10, 10, 10, 10, 10, 10)
prediction <- c(0.9, 0.8, 0.7, 0.5, 0.4, 0.3, 0.2)
NagelkerkeScaledR2(x, size, prediction)

# Using the example in fmsb::NagelkerkeR2
res <- glm(cbind(ncases,ncontrols) ~ agegp+alcgp+tobgp, data=esoph, family=binomial())
NagelkerkeScaledR2(x=esoph$ncases, size = esoph$ncases+esoph$ncontrols, 
                   prediction = res$fitted.values)
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
