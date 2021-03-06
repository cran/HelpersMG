% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScalePreviousPlot.R
\name{ScalePreviousPlot}
\alias{ScalePreviousPlot}
\title{Return the scale of the previous plot}
\usage{
ScalePreviousPlot()
}
\value{
A list with xlim and ylim
}
\description{
Return a list with the limits of the previous plot, the center, the range, and the position of label on this axe.
}
\details{
ScalePreviousPlot returns the scale of the previous plot
}
\examples{
\dontrun{
par(xaxs="i", yaxs="i")
plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
xlim= ScalePreviousPlot()$xlim[1:2]
ylim= ScalePreviousPlot()$ylim[1:2]
par(xaxs="r", yaxs="i")
plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
xlim= ScalePreviousPlot()$xlim[1:2]
ylim= ScalePreviousPlot()$ylim[1:2]
# Here is an example of the use of the label output
plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="", ylab="")
text(x=ScalePreviousPlot()$xlim["label"], y=ScalePreviousPlot()$ylim["center"], 
  xpd=TRUE, "Legend for Y axes", pos=3, srt=90)
text(x=ScalePreviousPlot()$xlim["center"], y=ScalePreviousPlot()$ylim["label"], 
  xpd=TRUE, "Legend for X axes", pos=1)
}
}
\seealso{
Other plot and barplot functions: 
\code{\link{barplot_errbar}()},
\code{\link{plot_add}()},
\code{\link{plot_errbar}()},
\code{\link{show_name}()}
}
\author{
Marc Girondot
}
\concept{plot and barplot functions}
