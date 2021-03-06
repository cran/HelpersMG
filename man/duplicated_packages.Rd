% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duplicated_packages.R
\name{duplicated_packages}
\alias{duplicated_packages}
\title{List the duplicated packages with their locations}
\usage{
duplicated_packages()
}
\value{
A data.frame with 4 elements for each duplicated packages:\cr
- versions: the version of the packages\cr
- libraries: the locations\cr
}
\description{
A data.frame with the duplicated packages and their locations and version.\cr
The columns Lib1 and Version1 should have the oldest version of the packages. Then 
you can try:\cr
li <- duplicated_packages()\cr
if (nrow(li) != 0) {\cr
    for (i in 1:nrow(li))\cr
        remove.packages(rownames(li)[i], lib=li[i, "Lib1"])\cr
}
}
\details{
duplicated_packages lists the duplicated packages with their locations
}
\examples{
\dontrun{
library(HelpersMG)
duplicated_packages()
}
}
\author{
Marc Girondot
}
