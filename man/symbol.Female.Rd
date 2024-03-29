% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.Female.R
\name{symbol.Female}
\alias{symbol.Female}
\title{Plot a female symbol in the plotting region}
\usage{
symbol.Female(centerx, centery, rayonx, lwd = 2, col = "black")
}
\arguments{
\item{centerx}{The x position of the center of the circle}

\item{centery}{The y position of the center of the circle}

\item{rayonx}{The size of the rayon in the scale of the x axis}

\item{lwd}{The width of the line of the symbol}

\item{col}{The color of the symbol}
}
\value{
Nothing
}
\description{
Plot a female symbol in the plotting region.
}
\details{
symbol.Female plot a female symbol in the plotting region
}
\examples{
\dontrun{
plot(x=1:2, y=c(10,20), type="n", bty="n", xlab="", ylab="")

rayonx <- 0.01
centerx <- 1.2
centery <- 15

symbol.Male(centerx=centerx, centery = centery, rayonx=rayonx)
symbol.Female(centerx=centerx+0.5, centery = centery, rayonx=rayonx)

rayonx <- 0.03
centerx <- 1.2
centery <- 18

symbol.Male(centerx=centerx, centery = centery, rayonx=rayonx, lwd=3)
symbol.Female(centerx=centerx+0.5, centery = centery, rayonx=rayonx, lwd=3, col="red")

rayonx <- 0.05
centerx <- 1.4
centery <- 13

symbol.Male(centerx=centerx, centery = centery, rayonx=rayonx, lwd=4, col="blue")
symbol.Female(centerx=centerx+0.5, centery = centery, rayonx=rayonx, lwd=4, col="red")
}
}
\seealso{
Other Symbol: 
\code{\link{symbol.Male}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{Symbol}
