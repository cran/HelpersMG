% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_ncdf.R
\name{format_ncdf}
\alias{format_ncdf}
\title{Return an array with ncdf data}
\usage{
format_ncdf(
  ncdf,
  label.latitude = "latitude",
  label.longitude = "longitude",
  label.time = "time",
  varid = NULL,
  longitude1 = NA,
  latitude1 = NA,
  longitude2 = NA,
  latitude2 = NA,
  package = "ncdf4",
  bathy = TRUE
)
}
\arguments{
\item{ncdf}{An object read from package ncdf4 or a file name of ncdf file}

\item{label.latitude}{Label of latitude}

\item{label.longitude}{Label of longitude}

\item{label.time}{Label of time}

\item{varid}{Name of variable to extract}

\item{longitude1}{Longitude for first corner}

\item{latitude1}{latitude for first corner}

\item{longitude2}{Longitude for second corner}

\item{latitude2}{latitude for second corner}

\item{package}{If ncdf is a file, give the package to use to open the file}

\item{bathy}{If TRUE, return a bathy object}
}
\value{
A list with two element: data is an array and time is the POSIX.lt time
}
\description{
Return a list with two elements: data is an array and time is the POSIX.lt time.\cr
Or if label.time is NULL or if bathy is TRUE, a bathy object.\cr
If varid is NULL, it shows the available variable and dimensions of the file.\cr
Bathymetry data can be download here: \cr
https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global
}
\details{
format_ncdf is used extract information from ncdf file
}
\examples{
\dontrun{
url <- "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/"
url <- paste0(url, "sst.day.mean.2012.v2.nc")
dest <- paste(Sys.getenv("HOME"), "/sst.day.mean.2012.v2.nc", sep="")
download.file(url, dest)
format_ncdf(dest)
}
}
\seealso{
Other ncdf: 
\code{\link{ind_long_lat}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{ncdf}
