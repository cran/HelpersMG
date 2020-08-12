#' format_ncdf is used extract information from ncdf file
#' @title Return an array with ncdf data
#' @author Marc Girondot
#' @return A list with two element: data is an array and time is the POSIX.lt time
#' @param ncdf An object read from package ncdf4 or a file name of ncdf file
#' @param label.latitude Label of latitude
#' @param label.longitude	Label of longitude
#' @param label.time	Label of time
#' @param varid	Name of variable to extract
#' @description Return a list with two element: data is an array and time is the POSIX.lt time./cr
#' If varid is NULL, it shows the available variable of the file.
#' @family ncdf
#' @examples
#' \dontrun{
#' url <- "ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/"
#' url <- paste0(url, "sst.day.mean.2012.v2.nc")
#' dest <- paste(Sys.getenv("HOME"), "/sst.day.mean.2012.v2.nc", sep="")
#' download.file(url, dest)
#' format_ncdf(dest)
#' }
#' @export


format_ncdf <- function(ncdf, 
                        label.latitude="latitude", 
                        label.longitude="longitude",
                        label.time="time", 
                        varid=NULL) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("ncdf4 package is necessary for this function")
  }
  if (!requireNamespace("RNetCDF", quietly = TRUE)) {
    stop("RNetCDF package is necessary for this function")
  }
  # library("ncdf4")
  if (class(ncdf) != "ncdf4")
    dta <- getFromNamespace("nc_open", ns="ncdf4")(ncdf)
  else
    dta <- ncdf
  
  if (is.null(varid)) {
    if (length(names(dta$var)) > 1)
      message(paste0("The available variables are: ", paste(names(dta$var), collapse = ", ")))
    else
      message(paste0("The available variable is: ", paste(names(dta$var), collapse = ", ")))
    return(invisible())
  }
  
  if (all(names(dta$dim) != label.latitude)) stop(paste0("Check the latitude label: ", paste(names(dta$dim), collapse = ", ")))
  if (all(names(dta$dim) != label.longitude)) stop(paste0("Check the longitude label: ", paste(names(dta$dim), collapse = ", ")))
  if (all(names(dta$dim) != label.time)) stop(paste0("Check the time label: ", paste(names(dta$dim), collapse = ", ")))
  
  carte3D <- getFromNamespace("ncvar_get", ns="ncdf4")(dta, varid=varid, 
                                                       start=c(1,1,1), 
                                                       count=c(dta$dim[[label.longitude]]$len, 
                                                               dta$dim[[label.latitude]]$len, 
                                                               dta$dim[[label.time]]$len))
  
  c3d <- array(data=carte3D[], dim=c(dta$dim[[label.longitude]]$len, 
                                     dta$dim[[label.latitude]]$len, 
                                     dta$dim[[label.time]]$len))
  
  # library(RNetCDF)
  date.char <- getFromNamespace("utcal.nc", ns="RNetCDF")(dta$dim[[label.time]]$units, dta$dim[[label.time]]$vals, type="s")
  date.POSIXlt <- strptime(date.char, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  dnames <- list(dta$dim[[label.longitude]]$vals, dta$dim[[label.latitude]]$vals, date.char)
  names(dnames) <- names(dta$dim)
  dimnames(c3d) <- dnames
  
  if (dta$dim[[label.time]]$len > 1) {
    if (date.POSIXlt[1] > date.POSIXlt[2]) {
      c3d <- c3d[, , dim(c3d)[3]:1, drop = FALSE]
    }
  }
  
  if (dta$dim[[label.latitude]]$len > 1) {
    if (as.numeric(dta$dim[[label.latitude]]$vals[1]) > dta$dim[[label.latitude]]$vals[2]) {
      c3d <- c3d[, dim(c3d)[2]:1, , drop = FALSE]
    }
  }
  
  if (dta$dim[[label.longitude]]$len > 1) {
    if (as.numeric(dta$dim[[label.longitude]]$vals[1]) > dta$dim[[label.longitude]]$vals[2]) {
      c3d <- c3d[dim(c3d)[1]:1, , , drop = FALSE]
    }
  }
  
  return(list(data=c3d, time=date.POSIXlt))
}