% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_AIC.R
\name{compare_AIC}
\alias{compare_AIC}
\title{Compares the AIC of several outputs}
\usage{
compare_AIC(
  ...,
  factor.value = 1,
  silent = FALSE,
  FUN = function(x) specify_decimal(x, decimals = 2)
)
}
\arguments{
\item{...}{Successive results to be compared as lists.}

\item{factor.value}{The $value of the list object is multiplied by factor.value to calculate AIC.}

\item{silent}{If TRUE, nothing is displayed.}

\item{FUN}{Function used to show values}
}
\value{
A list with DeltaAIC and Akaike weight for the models.
}
\description{
This function is used to compare the AIC of several outputs obtained with the same data but with different set of parameters.\cr
The parameters must be lists with $aic or $AIC or $value and $par elements or if AIC(element) is defined.\cr
if \code{$value} and \code{$par} are present in the object, the AIC is calculated as \code{2*factor.value*value+2*length(par)}. If \code{$value} is -log(likeihood), then factor.value must be 1 and if \code{$value} is log(likeihood), then factor.value must be -1.\cr
If several objects are within the same list, their AIC are summed.\cr
For example, compare_AIC(g1=list(group), g2=list(separe1, separe2)) can be used to compare a single model onto two different sets of data against each set of data fited with its own set of parameters.\cr
Take a look at \code{ICtab} in package \code{bbmle} which is similar.
}
\details{
compare_AIC compares the AIC of several outputs obtained with the same data.
}
\examples{
\dontrun{
library("HelpersMG")
# Here two different models are fitted
x <- 1:30
y <- rnorm(30, 10, 2)+log(x)
plot(x, y)
d <- data.frame(x=x, y=y)
m1 <- lm(y ~ x, data=d)
m2 <- lm(y ~ log(x), data=d)
compare_AIC(linear=m1, log=m2)
# Here test if two datasets can be modeled with a single model
x2 <- 1:30
y2 <- rnorm(30, 15, 2)+log(x2)
plot(x, y, ylim=c(5, 25))
plot_add(x2, y2, col="red")
d2 <- data.frame(x=x2, y=y2)
m1_2 <- lm(y ~ x, data=d2)
x_grouped <- c(x, x2)
y_grouped <- c(y, y2)
d_grouped <- data.frame(x=x_grouped, y=y_grouped)
m1_grouped <- lm(y ~ x, data=d_grouped)
compare_AIC(separate=list(m1, m1_2), grouped=m1_grouped)
}
}
\seealso{
Other AIC: 
\code{\link{ExtractAIC.glm}()},
\code{\link{FormatCompareAIC}()},
\code{\link{compare_AICc}()},
\code{\link{compare_BIC}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{AIC}
