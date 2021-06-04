#' tide.info gets the annual tide calendar for one particular location.
#' @title Annual tide calendar for one particular location
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Return a data.frame with annual tide calendar.
#' @param location Textual information
#' @param year Year to get the calendar
#' @param longitude Longitude to search for
#' @param latitude Latitude to search for
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
#' points(tide[1, "Longitude"], tide[1, "Latitude"], col="red")
#' points(-32, -4, col="blue")
#' axis(1)
#' axis(2, las=1)
#'      
#' # Show the locations with data    
#' library(maps)
#' map()
#' axis(1)
#' axis(2, las=1)
#' points(tide_location[, c("longitude")], tide_location[, c("latitude")], 
#'        pch=".", col="red", cex=2)
#' }
#' @export

tide.info <- function(location=NULL, year=2021, longitude=NULL, latitude=NULL) {
  
  tide_location <- tide_location
  
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
    
  } else {
    indLl <- sqrt((tide_location$latitude-latitude)^2+(tide_location$longitude-longitude)^2)
    pos <- which.min(abs(indLl))[1]
  }
  
  location <- tide_location$name[pos]
  message(paste0("The location ", location, " has been chosen."))
  print(tide_location[pos, c("name", "longitude", "latitude", "timezone", "country")])
  
  Begin <- paste0(as.character(year), "-01-01 00:00")
  End <- paste0(as.character(year), "-12-31 23:59")
  
  out <- system(command = paste0('curl --silent --data-urlencode "location=', location,'" --data-urlencode "begin=', Begin,'" --data-urlencode "end=', End,'" "http://134.158.74.46:20000/tide"'), intern = TRUE)
  
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
  
  out <- cbind(out, DateTime.UTC=convert.tz(out$DateTime.local, tz="UTC"))
  out <- cbind(out, Longitude=tide_location$longitude[pos])
  out <- cbind(out, Latitude=tide_location$latitude[pos])
  
  
  out <- cbind(out, Tide.meter=as.numeric(gsub(" ft", "", out$Tide))*0.3084)[, c("Location", "Longitude", "Latitude", 
                                                                                 "Phase", "DateTime.local", "DateTime.UTC", "Tide.meter")]
  out <- na.omit(out)
  
  return(out)
}
