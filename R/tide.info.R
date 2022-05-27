#' tide.info gets the annual tide calendar for one particular location.
#' @title Annual tide calendar for one particular location
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a data.frame with annual tide calendar.
#' @param location Textual information about location name
#' @param year Year to get the calendar
#' @param longitude Longitude to search for
#' @param latitude Latitude to search for
#' @param force.tide.height If FALSE, can return a current speed rather than tide height
#' @family Periodic patterns of indices
#' @description Annual tide information.\cr
#' The columns are: Location, Longitude, Latitude, Phase, DateTime.local, DateTime.UTC, Tide.meter\cr
#' This function uses an API linking xtide software 
#' (https://flaterco.com/xtide/) with tide.info() function.\cr
#' You must have a working internet connection for this function.
#' @keywords Tide
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' Location <- "Les Hattes"
#' Year <- 2010
#' tide <- tide.info(Location, Year)
#' plot(tide[, "DateTime.local"], tide[, "Tide.meter"], 
#'      type="l", bty="n", las=1, 
#'      main=tide[1, "Location"], 
#'      xlab=as.character(Year), ylab="Tide level in meter")
#' 
#' Location <- "Hawaii"
#' Year <- 2010
#' tide <- tide.info(Location, Year)
#' 
#' Location <- "Hanamaulu Bay, Kauai Island, Hawaii"
#' Year <- 2010
#' tide <- tide.info(Location, Year)
#' plot(tide[, "DateTime.local"], tide[, "Tide.meter"], 
#'      type="l", bty="n", las=1, 
#'      main=tide[1, "Location"], 
#'      xlab=as.character(Year), ylab="Tide level in meter")
#'      
#' tide <- tide.info(year=2010, longitude=-32, latitude=-4)
#' library(maps)
#' map(database = "world", regions = "Brazil", asp=1, 
#'     xlim=c(-80, -30), ylim=c(-33, 5))
#' points(tide[1, "Longitude"], tide[1, "Latitude"], col="red", pch=19)
#' points(-32, -4, col="blue", pch=19)
#' axis(1)
#' axis(2, las=1)
#' 
#' # Show the locations with data    
#' library(maps)
#' map(xlim=c(-180, 180), ylim=c(-90, 90))
#' title("Locations with harmonics data")
#' axis(1, at=seq(from=-180, to=180, by=45))
#' axis(2, las=1, at=seq(from=-90, to=90, by=15))
#' points(getFromNamespace(x="tide_location", ns="HelpersMG")[, c("longitude")], 
#'        getFromNamespace(x="tide_location", ns="HelpersMG")[, c("latitude")], 
#'        pch=".", col="red", cex=2)
#' # Another example
#' tikei_lon  <- (-144.5465183)
#' tikei_lat <- -14.9505897
#' Year <- 2021
#' tikei_tide <- tide.info(year=Year, longitude=tikei_lon, latitude=tikei_lat)
#' plot(tikei_tide[, "DateTime.local"], tikei_tide[, "Tide.meter"], 
#'      type="l", bty="n", las=1, 
#'      main=tikei_tide[1, "Location"], 
#'      xlab=as.character(Year), ylab="Tide level in meter")
#' ## Another one
#' tikei_lon <- (-75.56861111)
#' tikei_lat <- 39.50083333
#' Year <- 2012
#' tikei_tide <- tide.info(year=Year, longitude=tikei_lon, latitude=tikei_lat)
#' 
#' library(mapdata)
#' map('worldHires', xlim=c(-77, -74), ylim=c(37, 40))
#' points(x=tikei_lon, y=tikei_lat, pch=19, col="red", cex=1)
#' points(x=tikei_tide$Longitude[1], y=tikei_tide$Latitude[2], 
#'        pch=19, col="blue", cex=1)
#' 
#' par(mar=c(4, 4, 2, 2))
#' plot(tikei_tide$DateTime.local, tikei_tide$Tide.meter, type="l")
#' }
#' @export

tide.info <- function(location=NULL, 
                      year=2021, 
                      longitude=NULL, latitude=NULL, 
                      force.tide.height=TRUE) {
  
  tide_location <- getFromNamespace(x="tide_location", ns="HelpersMG")
  
  if (is.null(location) & (is.null(longitude) | is.null(latitude))) {
    stop("Location or longitude/latitude must be provided.")
  }
  
  if (!is.null(location)) {
    
    pos <- which(grepl(location, tide_location$index))
    
    if (identical(pos, integer(0))) {
      stop("No location is found")
    }
    
    if (length(pos) != 1) {
      if (all(tide_location$name[pos[1]] == tide_location$name[pos[-1]])) {
        pos_time <- which.max(tide_location$date_imported[pos])
        pos <- pos[pos_time[1]]
      } else {
        print(tide_location$name[pos])
        stop("More that one location are found; be more precise.")
      }
    }
    dist <- NULL
    
  } else {
    indLl <- sqrt((tide_location$latitude-latitude)^2+(tide_location$longitude-longitude)^2)
    o_pos <- order(indLl)
    pos <- 1
    if (force.tide.height) {
      while (grepl("urrent", tide_location$name[o_pos[pos]])) {
        pos <- pos + 1
      }
    }
    pos <- o_pos[pos]
    toRad <- pi/180
    p1 <- matrix(c(latitude=latitude, longitude=longitude), ncol = 2)* toRad
    p2 <- as.matrix(tide_location[pos, c("latitude", "longitude")])* toRad
    r <- 6378137
      p = cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2], as.vector(r))
      dLat <- p[, 4] - p[, 2]
      dLon <- p[, 3] - p[, 1]
      a <- (sin(dLat/2))^2 + cos(p[, 2]) * cos(p[, 4]) * (sin(dLon/2))^2
      a <- pmin(a, 1)
      dist <- as.vector(2 * atan2(sqrt(a), sqrt(1 - a)) * p[, 5])
  }
  
  location <- tide_location$name[pos]
  message(paste0("The location ", location, " has been chosen."))
  if (!is.null(dist)) message(paste0("The distance between your setup and ", location, " is ", dist, " meters (Haversine distance)."))
  print(tide_location[pos, c("name", "longitude", "latitude", "timezone", "country")])
  
  # tz <- tide_location[pos, "timezone"]
  
  Begin <- paste0(as.character(year), "-01-01 00:00")
  End <- paste0(as.character(year), "-12-31 23:59")
  
  file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".txt")
  
  com <- paste0('curl --silent -o ', file, ' --data-urlencode "location=', location,'" --data-urlencode "begin=', Begin,'" --data-urlencode "end=', End,'" "http://134.158.74.46:20000/tide"')
  out <- system(command = com, intern = FALSE)
  out <- readChar(file, nchars=file.info(file)$size)
  
  out <- strsplit(out, split = "\",\"")[[1]]
  out <- gsub("\\[\"", "", out)
  out <- gsub("\"\\]", "", out)
  out <- strsplit(out, split = ",")
  
  out_df <- data.frame(Location=character(), 
                       Date=character(), Time=character(), 
                       Tide=character(), Phase=character())
  
  for (i in 1:length(out)) out_df <- rbind(out_df, as.data.frame(matrix(out[[i]], nrow=1)))
  out <- out_df
  rm(out_df)
  colnames(out) <- c("Location", "Date", "Time", "Tide", "Phase")
  
  out <- cbind(out, DateTime=paste(out$Date, out$Time))
  
  loc <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "C") 
  
  out <- cbind(out, DateTime.local=strptime(out$DateTime, format = "%Y-%m-%d %I:%M %p", tz=tide_location$timezone[pos]))
  Sys.setlocale("LC_TIME", loc) 
  
  # if (is.null(tz)) {
    out <- cbind(out, DateTime.UTC=convert.tz(out$DateTime.local, tz="UTC"))
  # } else {
  #   out <- cbind(out, DateTime.UTC=convert.tz(out$DateTime.local, tz=tz))
  # }
  
  out <- cbind(out, Longitude=tide_location$longitude[pos])
  out <- cbind(out, Latitude=tide_location$latitude[pos])
  
  unit <- NULL
  td <- NULL
  if (any(grepl("ft", out$Tide))) {
    td <- as.numeric(gsub(" ft", "", out$Tide))*0.3084
    unit <- "Tide.meter"
  }
  if (any(grepl("m", out$Tide[1]))) {
    td <- as.numeric(gsub(" m", "", out$Tide))
    unit <- "Tide.meter"
  }
  if (any(grepl("kt", out$Tide[1]))) {
    td <- as.numeric(gsub(" kt", "", out$Tide))
    unit <- "Tide.knot"
  }
  
  td <- data.frame(td=td)
  colnames(td) <- unit
  out <- cbind(out, td)[, c("Location", "Longitude", "Latitude", "Phase", "DateTime.local", "DateTime.UTC", unit)]
  out <- na.omit(out)
  
  return(out)
}
