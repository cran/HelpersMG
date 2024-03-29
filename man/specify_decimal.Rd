% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specify_decimals.R
\name{specify_decimal}
\alias{specify_decimal}
\title{Return a number as character with specified number of decimals}
\usage{
specify_decimal(x, decimals = NULL, decimal.point = ".")
}
\arguments{
\item{x}{The numbers to be formated}

\item{decimals}{Number of decimals to print}

\item{decimal.point}{Character to be used as decimal point}
}
\value{
A character
}
\description{
Return a number as character with specified number of decimals.
If a is a matrix, it will return a matrix of the same size and the same attributes.
}
\details{
specify_decimals format a number with specified number of decimals
}
\examples{
specify_decimal(x=pi, decimal.point=".")
specify_decimal(x=pi, decimals=4, decimal.point=".")
specify_decimal(x=c(pi, exp(1)), decimals=3, decimal.point=",")
specify_decimal(x=c(pi, exp(1)), decimal.point=",")
specify_decimal(x=c(pi*10, pi, pi/10, pi/100, pi/1000))
specify_decimal(x=c(pi=pi), decimal.point=".")
specify_decimal(x=matrix(pi*1:4, ncol=2), decimal.point=".")
m <- matrix(pi*1:4, ncol=2)
rownames(m) <- c("A", "B")
colnames(m) <- c("C", "D")
specify_decimal(x=m, decimal.point=".")
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
