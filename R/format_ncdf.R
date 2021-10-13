#' format_ncdf is used extract information from ncdf file
#' @title Return an array with ncdf data
#' @author Marc Girondot
#' @return A list with two element: data is an array and time is the POSIX.lt time
#' @param ncdf An object read from package ncdf4 or a file name of ncdf file
#' @param label.latitude Label of latitude
#' @param label.longitude	Label of longitude
#' @param label.time	Label of time
#' @param varid	Name of variable to extract
#' @param longitude1 Longitude for first corner
#' @param latitude1 latitude for first corner
#' @param longitude2 Longitude for second corner
#' @param latitude2 latitude for second corner
#' @param package If ncdf is a file, give the package to use to open the file
#' @param bathy If TRUE, return a bathy object
#' @description Return a list with two elements: data is an array and time is the POSIX.lt time.\cr
#' Or if label.time is NULL or if bathy is TRUE, a bathy object.\cr
#' If varid is NULL, it shows the available variable and dimensions of the file.\cr
#' Bathymetry data can be download here: \cr
#' https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global
#' @family ncdf
#' @examples
#' \dontrun{
#' url <- "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/"
#' url <- paste0(url, "sst.day.mean.2012.v2.nc")
#' dest <- paste(Sys.getenv("HOME"), "/sst.day.mean.2012.v2.nc", sep="")
#' download.file(url, dest)
#' format_ncdf(dest)
#' }
#' @export


format_ncdf <- function(ncdf                         , 
                        label.latitude="latitude"    , 
                        label.longitude="longitude"  ,
                        label.time="time"            , 
                        varid=NULL                   , 
                        longitude1=NA                , 
                        latitude1=NA                 ,  
                        longitude2=NA                , 
                        latitude2=NA                 , 
                        package = "ncdf4"            , 
                        bathy = TRUE                 ) {
  
  package <- match.arg(arg=package, choices=c("ncdf4", "RNetCDF"))
  
  if (
    (
      (
        any(class(ncdf) == "ncdf4") 
      ) | 
      ( 
        (
          any(class(ncdf) == "character")
        ) & (
          package == "ncdf4"
        ) 
      ) 
    ) & 
    (!requireNamespace("ncdf4", quietly = TRUE) ) ) {
    
    stop("ncdf4 package is necessary for this function")
  }
  if (
    (
      (
        any(class(ncdf) == "NetCDF") 
      ) | 
      ( 
        (
          any(class(ncdf) == "character")
        ) & (
          package == "RNetCDF"
        ) 
      ) 
    ) & (!requireNamespace("RNetCDF", quietly = TRUE) ) ) {
    
    stop("RNetCDF package is necessary for this function")
  }
  
  # library("ncdf4")
  if (any(class(ncdf) == "character")) {
    if (package == "ncdf4") {
      ncdf <- getFromNamespace("nc_open", ns="ncdf4")(ncdf)
    } else {
      ncdf <- getFromNamespace("open.nc", ns="RNetCDF")(ncdf)
    }
  }
  
  if (any(class(ncdf) == "ncdf4")) {
    vr <- names(ncdf$var)
    dm <- names(ncdf$dim)
  } else {
    # Je suppose que je suis en RNetCDF
    info <- getFromNamespace("file.inq.nc", ns="RNetCDF")(ncdf)
    vr <- NULL
    for (i in (info$ndims):(info$nvars-1)) {
      vr <- c(vr, getFromNamespace("var.inq.nc", ns="RNetCDF")(ncdf, variable = i)$name)
    }
    dm <- NULL
    for (i in 0:(info$ndims-1)) {
      dm <- c(dm, getFromNamespace("dim.inq.nc", ns="RNetCDF")(ncdf, dimension = i)$name)
    }
  }
  
  message(paste0("The available variables are: ", paste(vr, collapse = ", ")))
  message(paste0("The available dimensions are: ", paste(dm, collapse = ", ")))
  
  if (is.null(varid)) {
    return(invisible())
  }
  
  if (!is.null(label.latitude)) 
    if (all(dm != label.latitude)) stop(paste0("Check the latitude label: ", paste(dm, collapse = ", ")))
  if (!is.null(label.longitude)) 
    if (all(dm != label.longitude)) stop(paste0("Check the longitude label: ", paste(dm, collapse = ", ")))
  if (!is.null(label.time)) 
    if (all(dm != label.time)) stop(paste0("Check the time label: ", paste(dm, collapse = ", ")))
  
  if (all(c(!is.na(longitude1), !is.na(latitude1), 
            !is.na(longitude1), !is.na(latitude1)))) {
    ind1 <- ind_long_lat(ncdf=ncdf, 
                         long = longitude1, 
                         lat = latitude1, 
                         label.longitude = label.longitude, 
                         label.latitude = label.latitude)
    ind2 <- ind_long_lat(ncdf=ncdf, 
                         long = longitude2, 
                         lat = latitude2, 
                         label.longitude = label.longitude, 
                         label.latitude = label.latitude)
    start_lg <- min(ind1["indice.long"], ind2["indice.long"])
    count_lg <- abs(ind2["indice.long"] - ind1["indice.long"]) + 1
    start_lt <- min(ind1["indice.lat"], ind2["indice.lat"])
    count_lt <- abs(ind2["indice.lat"] - ind1["indice.lat"]) + 1
  } else {
    start_lg <- 1
    start_lt <- 1  
    if (any(class(ncdf) == "ncdf4")) {
      count_lg <- ncdf$dim[[label.longitude]]$len
      count_lt <- ncdf$dim[[label.latitude]]$len
    } else {
      count_lg <- getFromNamespace("dim.inq.nc", ns="RNetCDF")(ncdf, dimension  = label.longitude)$length
      count_lt <- getFromNamespace("dim.inq.nc", ns="RNetCDF")(ncdf, dimension  = label.latitude)$length
    }
  }
  
  
  if ((!is.null(label.time)) | (!bathy)) {
    
    
    if (any(class(ncdf) == "ncdf4")) {
      
      start <- NULL
      count <- NULL
      for (dm_ec in dm) {
        if (dm_ec == label.longitude) {
          start <- c(start, start_lg)
          count <- c(count, count_lg)
        } else if (dm_ec == label.latitude) {
          start <- c(start, start_lt)
          count <- c(count, count_lt)
        } else if (dm_ec == label.time) {
          start <- c(start, 1)
          count <- c(count, ncdf$dim[[label.time]]$len)
        } else {
          start <- c(start, 1)
          count <- c(count, 1)
        }
      }
      
      
      carte3D <- getFromNamespace("ncvar_get", ns="ncdf4")(ncdf, varid=varid,
                                                           start=start,
                                                           count=count)
      
      c3d <- array(data=carte3D[], dim=count[c(which(dm == label.longitude), 
                                               which(dm == label.latitude), 
                                               which(dm == label.time))])
      
      
      if (isTRUE(requireNamespace("RNetCDF", quietly = TRUE)) ) {
        date.char <- getFromNamespace("utcal.nc", ns="RNetCDF")(ncdf$dim[[label.time]]$units, 
                                                                ncdf$dim[[label.time]]$vals, 
                                                                type="s")
        date.POSIXlt <- strptime(date.char, format="%Y-%m-%d %H:%M:%S", tz="UTC")
      } else {
        warning("The date and time information required the package RNetCDF being installed")
        date.POSIXlt <- NULL
        date.char <- ncdf$dim[[label.time]]$vals
      }
      dnames <- list(ncdf$dim[[label.longitude]]$vals[start_lg:(start_lg+count_lg-1)], 
                     ncdf$dim[[label.latitude]]$vals[start_lt:(start_lt+count_lt-1)], 
                     date.char)
      
      names(dnames) <- c(label.longitude, label.latitude, label.time)
      
      if (ncdf$dim[[label.time]]$len > 1) {
        if (date.POSIXlt[1] > date.POSIXlt[2]) {
          # Je dois inverser les noms
          c3d <- c3d[, , dim(c3d)[3]:1, drop = FALSE]
          dnames[[3]] <- rev(dnames[[3]])
        }
      }
      
      if (ncdf$dim[[label.latitude]]$len > 1) {
        if (as.numeric(ncdf$dim[[label.latitude]]$vals[1]) > ncdf$dim[[label.latitude]]$vals[2]) {
          # Je dois inverser les noms
          c3d <- c3d[, dim(c3d)[2]:1, , drop = FALSE]
          dnames[[2]] <- rev(dnames[[2]])
        }
      }
      
      if (ncdf$dim[[label.longitude]]$len > 1) {
        if (as.numeric(ncdf$dim[[label.longitude]]$vals[1]) > ncdf$dim[[label.longitude]]$vals[2]) {
          # Je dois inverser les noms
          c3d <- c3d[dim(c3d)[1]:1, , , drop = FALSE]
          dnames[[1]] <- rev(dnames[[1]])
        }
      }
      
    } else {
      
      start <- NULL
      count <- NULL
      for (dm_ec in dm) {
        if (dm_ec == label.longitude) {
          start <- c(start, start_lg)
          count <- c(count, count_lg)
        } else if (dm_ec == label.latitude) {
          start <- c(start, start_lt)
          count <- c(count, count_lt)
        } else if (dm_ec == label.time) {
          start <- c(start, 1)
          count <- c(count, getFromNamespace("dim.inq.nc", ns="RNetCDF")(ncdf, dimension  = label.time)$length)
        } else {
          start <- c(start, 1)
          count <- c(count, 1)
        }
      }
      
      
      
      carte3D <- getFromNamespace("var.get.nc", ns="RNetCDF")(ncfile=ncdf, variable=varid, 
                                                              start=start, 
                                                              count=count)
      
      c3d <- array(data=carte3D[], dim=count[c(which(dm == label.longitude), 
                                               which(dm == label.latitude), 
                                               which(dm == label.time))])
      
      ulat <- getFromNamespace("var.get.nc", ns="RNetCDF")(ncfile=ncdf, variable = label.latitude, start = start_lt, count = count_lt) 
      ulon <- getFromNamespace("var.get.nc", ns="RNetCDF")(ncfile=ncdf, variable = label.longitude, start = start_lg, count = count_lg) 
      utime <- getFromNamespace("var.get.nc", ns="RNetCDF")(ncfile=ncdf, variable = label.time, start = 1, count = getFromNamespace("dim.inq.nc", ns="RNetCDF")(ncdf, dimension  = label.time)$length)
      
      date.char <- getFromNamespace("utcal.nc", ns="RNetCDF")(unitstring=getFromNamespace("att.get.nc", ns="RNetCDF")(ncdf, variable=label.time, attribute="units"), 
                                                              utime, type="s")
      date.POSIXlt <- strptime(date.char, format="%Y-%m-%d %H:%M:%S", tz="UTC")
      
      dnames <- list(ulon, 
                     ulat, 
                     date.char)
      names(dnames) <- c(label.longitude, label.latitude, label.time)
    }
    
    
    
    dimnames(c3d) <- dnames
    
    return(list(data=c3d, time=date.POSIXlt))
    
  } else {
    
    
    
    if (any(class(ncdf) == "ncdf4")) {
      
      carte2D <- getFromNamespace("ncvar_get", ns="ncdf4")(ncdf, varid=varid, 
                                                           start=c(start_lg, start_lt), 
                                                           count=c(count_lg, 
                                                                   count_lt))
      c2d <- array(data=carte2D[], dim=c(count_lg, 
                                         count_lt))
      
      dnames <- list(ncdf$dim[[label.longitude]]$vals[start_lg:(start_lg+count_lg-1)], 
                     ncdf$dim[[label.latitude]]$vals[start_lt:(start_lt+count_lt-1)] )
      names(dnames) <- names(ncdf$dim)
      
      if (ncdf$dim[[label.latitude]]$len > 1) {
        if (as.numeric(ncdf$dim[[label.latitude]]$vals[1]) > ncdf$dim[[label.latitude]]$vals[2]) {
          # Je dois inverser les noms
          c2d <- c2d[, dim(c2d)[2]:1, drop = FALSE]
          dnames[[2]] <- rev(dnames[[2]])
        }
      }
      
      if (ncdf$dim[[label.longitude]]$len > 1) {
        if (as.numeric(ncdf$dim[[label.longitude]]$vals[1]) > ncdf$dim[[label.longitude]]$vals[2]) {
          # Je dois inverser les noms
          c2d <- c2d[dim(c2d)[1]:1, , drop = FALSE]
          dnames[[1]] <- rev(dnames[[1]])
        }
      }
      
    } else {
      
      carte2D <- getFromNamespace("var.get.nc", ns="RNetCDF")(ncfile=ncdf, variable=varid, 
                                                              start=c(start_lg, start_lt), 
                                                              count=c(count_lg, 
                                                                      count_lt))
      
      c2d <- array(data=carte2D[], dim=c(count_lg, 
                                         count_lt))
      
      ulon <- getFromNamespace("var.get.nc", ns="RNetCDF")(ncfile=ncdf, variable = label.longitude, start = start_lg, count = count_lg)
      ulat <- getFromNamespace("var.get.nc", ns="RNetCDF")(ncfile=ncdf, variable = label.latitude, start = start_lt, count = count_lt)
      
      dnames <- list(ulon, ulat)
      names(dnames) <- c(label.longitude, label.latitude)
    }
    
    dimnames(c2d) <- dnames
    
    if (bathy) class(c2d) <- "bathy"
    
    return(c2d)
  }
}
