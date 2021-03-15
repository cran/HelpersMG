#' ind_long_lat is used to manage ncdf information
#' @title Return or the index in ncdf object from lat/longitude or inverse
#' @author Marc Girondot
#' @return Or the index in ncdf object from lat/longitude or inverse
#' @param ncdf An object read from package ncdf4, ncdf or RNetCDF
#' @param long Longitude in decimal format
#' @param lat	Latitude in decimal format
#' @param indice.long	Index of longitude
#' @param indice.lat	Index of latitude
#' @param label.longitude Name of argument for longitude, default is lon
#' @param label.latitude Name of argument for latitude, default is lat
#' @description Return or the index in ncdf object from lat/longitude or reverse.
#' @family ncdf
#' @examples
#' \dontrun{
#' url <- "ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/"
#' url <- paste0(url, "sst.day.mean.2012.v2.nc")
#' dest <- paste(Sys.getenv("HOME"), "/sst.day.mean.2012.v2.nc", sep="")
#' download.file(url, dest)
#' library("ncdf4")
#' dta2012 <- nc_open(dest)
#' indices <- ind_long_lat(ncdf=dta2012, lat=5.89, long=-20.56)
#' coordinates <- ind_long_lat(ncdf=dta2012, indice.lat=20, indice.long=30)
#' # library("RNetCDF")
#' # dta2012 <- open.nc(dest)
#' # indices <- ind_long_lat(ncdf=dta2012, lat=5.89, long=-20.56)
#' # coordinates <- ind_long_lat(ncdf=dta2012, indice.lat=20, indice.long=30)
#' # ncdf library is depreciated in CRAN
#' # library("ncdf")
#' # dta2012 <- open.ncdf(dest)
#' # indices <- ind_long_lat(ncdf=dta2012, lat=5.89, long=-20.56)
#' # coordinates <- ind_long_lat(ncdf=dta2012, indice.lat=20, indice.long=30)
#' }
#' @export


ind_long_lat <- function (ncdf = stop("The ncdf data must be supplied") , 
                          long = NULL                                   , 
                          lat = NULL                                    , 
                          indice.long = NULL                            , 
                          indice.lat = NULL                             , 
                          label.longitude = "lon"                       , 
                          label.latitude = "lat"                        ) 
{
  maxindicelt <- NULL
  maxindicelg <- NULL
  # if (class(ncdf) == "ncdf4") {
  #   if (!requireNamespace("ncdf", quietly = TRUE)) {
  #     stop("ncdf package is necessary for this function")
  #   }
  #   maxindicelt <- ncdf$dim[[label.latitude]]$len
  #   lt <- ncdf$dim[[label.latitude]]$vals
  #   lg <- ncdf$dim[[label.longitude]]$vals
  #   maxlt <- lt[maxindicelt]
  #   minlt <- lt[1]
  #   maxindicelg <- ncdf$dim[[label.longitude]]$len
  #   maxlg <- lg[maxindicelg]
  #   minlg <- lg[1]
  # }
  if (any(class(ncdf) == "ncdf4")) {
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("ncdf4 package is necessary for this function")
    }
    maxindicelt <- ncdf$dim[[label.latitude]]$len
    lt <- ncdf$dim[[label.latitude]]$vals
    maxlt <- lt[maxindicelt]
    minlt <- lt[1]
    maxindicelg <- ncdf$dim[[label.longitude]]$len
    lg <- ncdf$dim[[label.longitude]]$vals
    maxlg <- lg[maxindicelg]
    minlg <- lg[1]
  }
  if (any(class(ncdf) == "NetCDF")) {
    if (!requireNamespace("RNetCDF", quietly = TRUE)) {
      stop("RNetCDF package is necessary for this function")
    }
    maxindicelt <- getFromNamespace("dim.inq.nc", ns = "RNetCDF")(ncfile = ncdf, 
                                                                  label.latitude)$length
    lt <- getFromNamespace("var.get.nc", ns = "RNetCDF")(ncfile = ncdf, variable = label.latitude)
    maxlt <- lt[maxindicelt]
    minlt <- lt[1]
    maxindicelg <- getFromNamespace("dim.inq.nc", ns = "RNetCDF")(ncfile = ncdf, 
                                                                  label.longitude)$length
    lg <- getFromNamespace("var.get.nc", ns = "RNetCDF")(ncfile = ncdf, 
                                                         variable = label.longitude)
    maxlg <- lg[maxindicelg]
    minlg <- lg[1]
  }
  if (is.null(maxindicelt) | is.null(maxindicelg)) {
    stop("Check the ncdf data; it is not recognized")
  }
  if (!is.null(long) & !is.null(lat)) {
    # which.min(abs(lt - lat))
    if ((long > maxlg) | (long < minlg)) long <- long%%360
    if ((lat > maxlt) | (lat < minlt)) lat <- ((lat + 90)%%180) - 90
    return(c(indice.long = which.min(abs(lg - long)), indice.lat = which.min(abs(lt - lat))))
  } else {
    if (!is.null(indice.long) & !is.null(indice.lat)) {
      return(c(long = lg[indice.long], lat = lt[indice.lat]))
    }
    else {
      stop("Check the parameters")
    }
  }
}
