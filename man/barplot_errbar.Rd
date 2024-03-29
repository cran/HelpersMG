% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barplot_errbar.R
\name{barplot_errbar}
\alias{barplot_errbar}
\title{Plot a barplot graph with error bar on y}
\usage{
barplot_errbar(
  ...,
  errbar.y = NULL,
  errbar.y.plus = NULL,
  errbar.y.minus = NULL,
  y.plus = NULL,
  y.minus = NULL,
  errbar.tick = 1/50,
  errbar.lwd = par("lwd"),
  errbar.lty = par("lty"),
  errbar.col = par("fg"),
  add = FALSE
)
}
\arguments{
\item{...}{Parameters for barplot() such as main= or ylim=}

\item{errbar.y}{The length of error bars for y. Recycled if necessary.}

\item{errbar.y.plus}{The length of positive error bars for y. Recycled if necessary.}

\item{errbar.y.minus}{The length of negative error bars for y. Recycled if necessary.}

\item{y.plus}{The absolut position of the positive error bar for y. Recycled if necessary.}

\item{y.minus}{The absolut position of the nagative error bar for y. Recycled if necessary.}

\item{errbar.tick}{Size of small ticks at the end of error bars defined as a proportion of total width or height graph size.}

\item{errbar.lwd}{Error bar line width, see par("lwd")}

\item{errbar.lty}{Error bar line type, see par("lwd")}

\item{errbar.col}{Error bar line color, see par("col")}

\item{add}{If true, add the graph to the previous one.}
}
\value{
A numeric vector (or matrix, when beside = TRUE), say mp, giving the coordinates of all the bar midpoints drawn, useful for adding to the graph.\cr
If beside is true, use colMeans(mp) for the midpoints of each group of bars, see example.
}
\description{
To plot data, just use it as a normal barplot but add the errbar.y
values or errbar.y.minus, errbar.y.plus if bars for y axis are
asymetric. Use y.plus and y.minus to set absolut limits for
error bars. Note that y.plus and y.minus have priority over errbar.y,
errbar.y.minus and errbar.y.plus.
}
\details{
barplot_errbar plot a barplot with error bar on y
}
\examples{
\dontrun{
barplot_errbar(rnorm(10, 10, 3), 
	xlab="axe x", ylab="axe y", bty="n", 
		errbar.y.plus=rnorm(10, 1, 0.1), col=rainbow(10), 
		names.arg=paste("Group",1:10), cex.names=0.6)
y <- rnorm(10, 10, 3)
barplot_errbar(y, 
               	xlab="axe x", ylab="axe y", bty="n", 
            		y.plus=y+2)
}
}
\seealso{
\code{plot_errorbar}

Other plot and barplot functions: 
\code{\link{ScalePreviousPlot}()},
\code{\link{plot_add}()},
\code{\link{plot_errbar}()},
\code{\link{show_name}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{plot and barplot functions}
